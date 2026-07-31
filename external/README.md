# External dependencies

Place the downloaded dependency folders in this directory, keeping these
folder names:

```text
external/
├─ slepian_alpha-master/
├─ slepian_bravo-master/
├─ slepian_delta-master/
├─ slepian_zero-master/
└─ GRACE-filter-master/    # optional; required for DDK products
```

The same parent directory can be stored elsewhere and passed as the second
argument to `startup_STPC_GA`, or selected with the `STPC_DEPENDENCIES`
environment variable.

Third-party source trees are intentionally not duplicated in this
repository. See the main README for download links and licenses.
