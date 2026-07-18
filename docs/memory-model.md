# Memory model

Two objects own all memory; nothing is allocated during event processing.

**`ilreco_config`** — the calorimeter: geometry, shower-profile tables,
optional cell mask, tunables. Created once in `ilreco_config_create`,
immutable after.

**`ilreco_workspace`** — the algorithm's working memory: hit/address
buffers, clustering scratch, per-peak weight tables, output staging. All of
it is **one memory block**, allocated in `ilreco_workspace_create` and
sized for the config's grid. Every event reuses that block:
`ilreco_reconstruct` allocates **nothing**, so per-event cost is flat
(~40 µs for a single shower).

## Ownership and destruction order

A workspace references its config; the config must outlive it. Free in
reverse order of creation:

```c
ilreco_workspace_destroy(ws);    /* every workspace first */
ilreco_config_destroy(cfg);      /* the config last       */
```

`destroy(NULL)` is a no-op. Destroying a config does not free its
workspaces: you free what you created.

## In C++

Scope enforces the order — declare the `Config` before the `Workspace` and
the destructors run in the correct reverse order:

```cpp
{
    ilreco::Config config(20, 20, "prof_pwo.dat");
    ilreco::Workspace workspace(config);
    auto clusters = workspace.reconstruct(hits);
}   // workspace destroyed first, then config — always correct
```

Both classes are move-only, so ownership can be transferred (into a vector
of per-thread workspaces, for example) but never duplicated.
