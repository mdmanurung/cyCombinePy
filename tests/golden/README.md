# R cyCombine Golden Fixtures

The R parity tests look for locally generated fixtures under versioned
directories:

```text
tests/golden/cycombine_r_<version>/
  manifest.json
  input.csv
  normalize_scale.npz
  normalize_rank.npz
  emd.csv
  mad.csv
  corrected_fixed_labels.npz
```

`manifest.json` must describe the exact environment and files used to create
the fixtures. It should include:

- `schema_version`
- `r_version`
- `cycombine.version`
- `cycombine.commit`
- `package_versions`
- `random_seed`
- `provenance`
- `sha256` hashes for every generated fixture file

These files are generated locally from an installed R `cyCombine` environment.
They are not downloaded from the network and should not depend on network
access during generation.
