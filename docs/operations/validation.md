# Validate and publish

Validation protects both scientific usefulness and the public/private
boundary.

## Completion checks

A release cannot be sealed when:

- a participant file is unexpectedly missing, zero bytes, malformed or has a
  checksum mismatch;
- manifests, sample sheet and answer key disagree about membership;
- a scored answer is missing or pending;
- an answer lies outside the submission schema;
- a source identity remains in a public sequence header;
- a faulted sample lacks passing materialisation evidence; or
- a required contract artifact is absent.

Run both checks after SLURM completes:

```bash
pixi run genomepuzzle release validate \
  --release-dir generated/round/release \
  --require-complete

pixi run genomepuzzle release inspect \
  --release-dir generated/round/release
```

`inspect` summarises the release type, sample count, participant file count,
normal/troublesome balance and complete-bundle SHA-256 digest.

## Scientific review

Automated validation is necessary but not sufficient. Before publication:

1. Review `private/implant_manifest.json`, including every `fault_type`.
2. Review `private/validation_report.json`.
3. Confirm expected answers were obtained from final participant files.
4. Inspect every troublesome sample and a representative clean sample.
5. Confirm the difficulty and pass threshold match the intended audience.
6. Confirm no public text, filename or sequence header exposes source identity.

## Publishing boundary

Hand the completed release directory to a contract-aware publisher. The
publisher should:

- verify `COMPLETE.json` and the complete-bundle digest;
- upload only `public/` to participant-facing storage;
- keep private truth in access-controlled server storage;
- import `release.json`, `submission_schema.json` and
  `scoring_policy.json`; and
- upload the public manifest last, so a partial upload is never advertised as
  ready.

Never publish a release with a plain recursive copy of its root directory.
