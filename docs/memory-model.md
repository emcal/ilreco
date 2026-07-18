# Memory model

Two objects own all memory; nothing is allocated during event processing.

**`ilreco_config`** — the calorimeter: geometry, shower-profile tables,
optional cell mask, tunables. Everything it will ever need is allocated in
`ilreco_config_create`; after that it is immutable and shared.

**`ilreco_workspace`** — the algorithm's working memory: hit/address
buffers, clustering scratch, per-peak weight tables, output staging. All of
it lives in **one allocation** (an internal arena), sized for the config's
grid at `ilreco_workspace_create`. Every event reuses the same arena —
`ilreco_reconstruct` performs **zero allocations**, which is what makes
per-event cost flat (~40 µs for a single shower) and the library usable in
tight online loops.

## Ownership and destruction order

A workspace references its config; the config must outlive it. Free in
reverse order of creation:

```c
ilreco_workspace_destroy(ws);    /* every workspace first */
ilreco_config_destroy(cfg);      /* the config last       */
```

`destroy(NULL)` is a no-op, and destroying a config never touches its
workspaces — there is no hidden registry; you free what you created.

## In C++

Scope does the ordering for you — declare the `Config` before the
`Workspace` and the destructors run in the correct reverse order:

```cpp
{
    ilreco::Config config(20, 20, "prof_pwo.dat");
    ilreco::Workspace workspace(config);
    auto clusters = workspace.reconstruct(hits);
}   // workspace destroyed first, then config — always correct
```

Both classes are move-only, so ownership can be transferred (into a vector
of per-thread workspaces, for example) but never duplicated.
