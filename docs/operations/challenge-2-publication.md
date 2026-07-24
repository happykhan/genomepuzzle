# Challenge 2 publication receipt

This is the durable publication record for the contract 2.1 practice and
Challenge 2 datasets uploaded on 24 July 2026.

## Source revisions

The final bundles accurately record two generator revisions:

- typing and hybrid: `433003951403756f4d786f3c4a81780ae574a0fd`;
- assembly and outbreak: `eb420a275a042dcf71c0f663fda35ad56eb1a83a`.

Assembly and outbreak were deliberately regenerated at the later revision to
use parameterised catastrophic read-pair counts. The annotated
`challenge-2-release` Git tag identifies the final tracked release record and
documents these per-exercise generator revisions; it does not pretend all
eight bundles were generated from one commit.

## Release digests

| Release ID | Bundle SHA-256 |
| --- | --- |
| `2026-website-typing-practice` | `c9fcdb57e1648e58d660a49d48399654da363e8c0d53f6ca2a5c4434ce020c7b` |
| `2026-website-assembly-practice` | `871188425151ec40cd88d3a45d5ea86f896802b7d1471bea3beca4324e5e371e` |
| `2026-website-hybrid-practice` | `464de50827d1e86eea482516eb574f499629560cb8e90369efc7df2d9efee0b9` |
| `2026-website-outbreak-practice` | `e61bc3ae65de55cd834a0b01fc523d40bef01a225ceda244750dadddd29952d6` |
| `challenge-2-typing` | `8307a36d4b979293ee2ea084b0a68278dc754ebcdd2a3d27fe34b7a38a7d6a12` |
| `challenge-2-assembly` | `6b9fef4a60d02fa4d4e5789002c38e2dfd37e85cde4c21037d62f490efca88e7` |
| `challenge-2-hybrid` | `add243400399180bfb31a2ea52b50b9df95bd98f7b72b516aaadd7b79a6bc506` |
| `challenge-2-outbreak` | `7311cc161a54ccba696716ddc1da0912480eee875eb8ad49bd80877ec989d0b2` |

The bundle digest covers the sealed local release contract. R2 additionally
stores a URL-enriched `dataset_manifest.json` for website consumption.

## R2 layout

Practice participant objects are in `ghrupuzzle-practice`. Challenge inputs
must be protected until the configured assessment window, so Challenge 2
participant objects are in `ghrupuzzle-private` alongside private truth.
Practice truth is also stored in the private bucket.

Every prefix uses:

```text
releases/{release_id}/{exercise}/{mode}
```

Post-publication inventory:

| Bucket | Objects | Bytes |
| --- | ---: | ---: |
| `ghrupuzzle-practice` | 69 | 1,735,909,811 |
| `ghrupuzzle-private` | 175 | 4,661,791,003 |

The private inventory consists of 122 Challenge 2 participant/contract
objects (4,661,554,925 bytes) and 53 private truth/provenance objects
(236,078 bytes). It contains no participant submissions or certificates.

## Replacement and verification

Desired contract 2.1 objects were uploaded or overwritten before stale-object
deletion. One obsolete practice object and four obsolete private objects were
then deleted.

Verification confirmed:

- exact desired key sets in both buckets;
- matching remote object sizes and SHA-256 metadata for all 244 objects;
- all eight remote manifests report schema `2.1`, the expected exercise, mode,
  release ID and sample count;
- Challenge 2 zero-byte assembly/read objects have remote size zero;
- deliberately missing hybrid long-read and outbreak R2 keys are absent;
- no challenge prefix or `/private/` key exists in the practice bucket.

Website registration remains a separate deployment action. GHRU Puzzles must
support contract 2.1 before these uploaded releases are registered.
