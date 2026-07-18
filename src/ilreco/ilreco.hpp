/*
 * ilreco.hpp — header-only C++ RAII wrapper over the ilreco C API.
 *
 * The C API (ilreco.h) stays the ground truth; this wrapper only manages
 * lifetimes and translates errors into exceptions:
 *
 *   ilreco::Config config(20, 20, "prof_pwo.dat");   // once per process
 *   ilreco::Workspace workspace(config);             // once per thread
 *   std::vector<ilreco_cluster> clusters =
 *       workspace.reconstruct({{12, 14, 1.85}, {13, 14, 0.42}});
 *
 * Threading rule (same as the C API): the Config is immutable and shared,
 * one Workspace per thread, and the Config must outlive every Workspace
 * created from it.
 */
#ifndef CALORIMETRY_STUDIES_ILRECO_HPP
#define CALORIMETRY_STUDIES_ILRECO_HPP

#include <ilreco.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace ilreco {

/* Immutable calorimeter configuration: geometry + shower-profile tables +
 * tunables. Create once, tune before the first Workspace, share across
 * threads. Move-only. */
class Config {
public:
    Config(std::int32_t n_cols, std::int32_t n_rows,
           const std::string& profile_path) {
        char error[256] = {};
        config_ = ilreco_config_create(n_cols, n_rows, profile_path.c_str(),
                                       error, sizeof error);
        if (!config_) {
            throw std::runtime_error(std::string("ilreco config: ") + error);
        }
    }
    ~Config() { ilreco_config_destroy(config_); }

    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
    Config(Config&& other) noexcept : config_(other.config_) {
        other.config_ = nullptr;
    }
    Config& operator=(Config&& other) noexcept {
        if (this != &other) {
            ilreco_config_destroy(config_);
            config_ = other.config_;
            other.config_ = nullptr;
        }
        return *this;
    }

    /* Cell-existence mask for non-rectangular shapes: row-major
     * n_rows*n_cols entries, nonzero = the cell physically exists.
     * Call before the first Workspace is created. */
    void set_cell_mask(const std::vector<unsigned char>& mask) {
        if (ilreco_config_set_cell_mask(config_, mask.data()) != 0) {
            throw std::invalid_argument("ilreco: invalid cell mask");
        }
    }
    void set_seed_threshold(double min_seed_gev) {
        ilreco_config_set_seed_threshold(config_, min_seed_gev);
    }
    void set_zcal(double zcal_cm) { ilreco_config_set_zcal(config_, zcal_cm); }
    void set_hole_classification(bool enabled) {
        ilreco_config_set_hole_classification(config_, enabled ? 1 : 0);
    }

    const ilreco_config* get() const { return config_; }

private:
    ilreco_config* config_ = nullptr;
};

/* Per-thread working memory. One Workspace per thread, never shared;
 * reused for every event (reconstruction allocates nothing). Move-only. */
class Workspace {
public:
    explicit Workspace(const Config& config) {
        workspace_ = ilreco_workspace_create(config.get());
        if (!workspace_) {
            throw std::runtime_error("ilreco: workspace allocation failed");
        }
    }
    ~Workspace() { ilreco_workspace_destroy(workspace_); }

    Workspace(const Workspace&) = delete;
    Workspace& operator=(const Workspace&) = delete;
    Workspace(Workspace&& other) noexcept : workspace_(other.workspace_) {
        other.workspace_ = nullptr;
    }
    Workspace& operator=(Workspace&& other) noexcept {
        if (this != &other) {
            ilreco_workspace_destroy(workspace_);
            workspace_ = other.workspace_;
            other.workspace_ = nullptr;
        }
        return *this;
    }

    /* Reconstruct one event; returns clusters energy-descending.
     * Throws std::invalid_argument on invalid input (hit outside the grid
     * or on a masked-out cell). max_clusters bounds the readout, not the
     * reconstruction. */
    std::vector<ilreco_cluster> reconstruct(const std::vector<ilreco_hit>& hits,
                                            std::int32_t max_clusters = 128) {
        std::vector<ilreco_cluster> clusters(max_clusters);
        const std::int32_t n_found = ilreco_reconstruct(
            workspace_, hits.data(), static_cast<std::int32_t>(hits.size()),
            clusters.data(), max_clusters);
        if (n_found < 0) {
            throw std::invalid_argument(
                "ilreco: invalid hit (outside the grid or on a masked-out cell)");
        }
        clusters.resize(n_found < max_clusters ? n_found : max_clusters);
        return clusters;
    }

    ilreco_workspace* get() { return workspace_; }

private:
    ilreco_workspace* workspace_ = nullptr;
};

}  // namespace ilreco

#endif /* CALORIMETRY_STUDIES_ILRECO_HPP */
