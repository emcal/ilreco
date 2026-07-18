# Advanced multithreading: when threads are not yours

The [threading model](./threading) is trivial when you own the worker
loop: create one workspace per thread. Frameworks often take that control
away — a thread pool, a task scheduler, or an event-processing framework
like JANA2 calls *your* code from *its* threads, and a factory has no idea
which thread invokes it, or how many exist.

Two patterns cover this. Both keep the guarantee: results are bitwise
independent of which thread ran which event.

## Pattern 1 — workspace pool (recommended)

A mutex-protected free-list of workspaces. Each call checks one out,
reconstructs, checks it back in. The pool grows to the peak concurrency and
no further; destruction is trivial because the pool owns everything.

```cpp
#include <ilreco.hpp>

#include <memory>
#include <mutex>
#include <vector>

// Reconstruction service callable from ANY thread (e.g. a JANA2 factory).
class ReconstructionService {
public:
    explicit ReconstructionService(int n_cols, int n_rows,
                                   const std::string& profile)
        : config_(n_cols, n_rows, profile) {}

    std::vector<ilreco_cluster> reconstruct(const std::vector<ilreco_hit>& hits) {
        std::unique_ptr<ilreco::Workspace> workspace = acquire();
        auto clusters = workspace->reconstruct(hits);
        release(std::move(workspace));
        return clusters;
    }

private:
    std::unique_ptr<ilreco::Workspace> acquire() {
        {
            std::lock_guard<std::mutex> lock(pool_mutex_);
            if (!pool_.empty()) {
                auto workspace = std::move(pool_.back());
                pool_.pop_back();
                return workspace;
            }
        }
        return std::make_unique<ilreco::Workspace>(config_);
    }
    void release(std::unique_ptr<ilreco::Workspace> workspace) {
        std::lock_guard<std::mutex> lock(pool_mutex_);
        pool_.push_back(std::move(workspace));
    }

    ilreco::Config config_;             // declared first: destroyed last
    std::mutex pool_mutex_;
    std::vector<std::unique_ptr<ilreco::Workspace>> pool_;
};
```

Cost: one mutex lock/unlock per event — nanoseconds against ~40 µs of
reconstruction. The member order makes destruction correct by construction:
the pool (workspaces) is destroyed before `config_`.

## Pattern 2 — `thread_local`

```cpp
std::vector<ilreco_cluster> reconstruct(const ilreco::Config& config,
                                        const std::vector<ilreco_hit>& hits) {
    thread_local ilreco::Workspace workspace(config);
    return workspace.reconstruct(hits);
}
```

Shorter, and lock-free — but with two caveats the pool doesn't have: the
workspace binds to the *first* config the thread sees (wrong with several
calorimeters in one program), and it is destroyed only at thread exit,
which on some runtimes is after your config is gone. Prefer the pool
unless the program provably has one calorimeter and framework-managed
thread lifetimes.

## This is exactly what the Python binding does

`ilreco.Calorimeter` embeds pattern 1 in its C core: any number of Python
threads — or the internal `n_jobs` chunk workers — may call
`reconstruct()` on the same object concurrently, each call borrowing a
workspace from the pool. See [Python API](./python-api).
