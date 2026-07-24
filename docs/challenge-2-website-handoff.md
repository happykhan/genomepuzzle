# Challenge 2 website handoff

This document is the implementation and acceptance-test brief for delivering
GenomePuzzle contract 2.1 releases through GHRU Puzzles.

## Release inventory

The website must consume these exact release IDs. The Challenge 2 IDs are not
year-labelled.

| Exercise | Practice release ID | Challenge release ID |
|---|---|---|
| Typing | `2026-website-typing-practice` | `challenge-2-typing` |
| Short-read assembly | `2026-website-assembly-practice` | `challenge-2-assembly` |
| Hybrid assembly | `2026-website-hybrid-practice` | `challenge-2-hybrid` |
| Outbreak investigation | `2026-website-outbreak-practice` | `challenge-2-outbreak` |

Every release uses schema version `2.1`. Practice public assets belong in the
practice R2 bucket. Challenge public assets and all private truth belong in the
private R2 bucket. Private files must never be returned by participant download
routes.

The object prefix is:

```text
releases/{release_id}/{exercise}/{mode}
```

## Why the website needs contract-driven UX

The datasets deliberately include missing, empty, truncated, contaminated, and
discordant inputs. These are assessment cases, not failed uploads. The website
must present exactly what the release manifest contains without trying to
"repair" it or rejecting a release because a conventional input role is absent.

Independent blinded analysis recovered all 48 QC decisions and all 48 failure
reasons exactly. The scientific fault model is therefore sufficiently clear.
The remaining risk is the website accidentally making valid answers difficult
to express, scoring unavailable fields, or revealing that a sample was
deliberately modified.

## Required release compatibility

1. Accept `schema_version: "2.1"` in release registration, retrieval,
   submission parsing, and scoring. Any hard-coded `2.0` check must be updated
   or replaced with an explicit supported-version set.
2. Load fields, labels, descriptions, allowed values, normalizers,
   `required_when`, scorer type, and score conditions from each release's
   `submission_schema.json` and private `scoring_policy.json`. Do not maintain a
   second exercise-specific field contract in the frontend.
3. Treat the manifest as authoritative for available files. A missing R1, R2,
   or long-read file can be intentional. A zero-byte object is also a valid
   downloadable assessment input.
4. Require exactly the release sample IDs in a submission. Reject duplicate,
   missing, and unexpected sample IDs with clear row-level messages.
5. Keep public manifests, schemas, instructions, and participant files
   separate from `private/answer_key.json`, implant manifests, provenance, and
   validation reports.

## Submission form

The first fields shown for every sample should be:

```text
sample_id, qc_status, failure_reason
```

`sample_id` is fixed by the release and must not be editable.

`qc_status` is a controlled choice containing only `PASS` and `FAIL`.

`failure_reason` is a controlled choice populated from the release schema:

- a passing row must use `NONE`;
- a failing row must use one permitted categorical failure reason;
- changing QC to `PASS` should automatically select `NONE`;
- changing QC to `FAIL` should require a non-`NONE` reason.

Do not use a free-text failure field. The purpose of the small vocabulary is to
make the decision objective and machine-scoreable.

### Conditional fields

When `qc_status` is `PASS`, show and require the scored analytical fields for
that exercise.

When `qc_status` is `FAIL`:

- require only `qc_status` and `failure_reason`;
- hide or visibly disable analytical fields that are unavailable;
- do not require species, typing, assembly, or cluster results;
- do not score those analytical fields, even if a participant entered a value.

This follows `score_when` in `scoring_policy.json`. It prevents participants
being penalised for refusing to analyse data they have correctly rejected.

Unscored fields such as assembler, assembly statistics, and notes may remain
available for evidence and manual review, but must be clearly labelled as
unscored.

### Exercise-specific fields

- Typing: species, ST, K locus, capsule type, wzi, O locus, O type, and
  carbapenemase.
- Short-read assembly: species is scored for passing samples; assembler and
  assembly statistics are supporting evidence.
- Hybrid assembly: species is scored for passing samples; assembler and
  assembly statistics are supporting evidence.
- Outbreak: species and cluster are scored for passing samples.

Do not add genotyping fields to assembly, hybrid, or outbreak submissions.
Do not add chromosome/plasmid status to the hybrid submission.

## Normalisation and scoring details

1. Apply the normalizer declared for each field before comparison.
2. Sequence type should accept equivalent human forms such as `14` and `ST14`.
3. Species and case-insensitive fields should tolerate harmless whitespace and
   case differences.
4. Carbapenemase values use unordered-list comparison; list order must not
   affect the score.
5. Kleborate uses the literal `-` for an unavailable `wzi`. A blinded analyst
   naturally submitted a blank instead, producing the only field mismatch in
   180 non-partition comparisons. Either:
   - normalise blank and `-` to the same unavailable value for fields where the
     answer key is `-`; or
   - use a controlled UI value such as `Not available (-)` that serialises to
     `-`.
6. Outbreak cluster labels are arbitrary. `cluster-a`, `1`, and `red` are
   equally valid labels if they describe the correct partition. Score whether
   each pair of passing samples is placed together or apart, not the literal
   cluster text.
7. A correctly failed outbreak sample must be excluded from partition scoring.
8. Calculate the final score from the private scoring policy and retain
   field-level scoring evidence for review.

## Participant feedback

Before submission, show:

- the expected columns;
- allowed controlled values;
- which fields are scored;
- which fields become required for `PASS`;
- a downloadable CSV template;
- validation errors attached to the affected row and field.

After a practice submission:

- score immediately;
- show per-sample QC and failure-reason feedback;
- explain analytical-field errors only for samples expected to pass;
- allow the participant to retry at any time;
- make the practice answer/explanation available according to the chosen
  practice feedback policy.

During a timed challenge:

- allow downloads and submissions only inside the configured challenge window;
- allow participants to replace their own submission before closing;
- preserve every submission version for audit;
- do not reveal private truth or detailed correctness before the challenge
  feedback time;
- show a clear receipt containing submission time and version.

## Manual review

Participants need a `Request review` action with an optional explanation.
Administrators need a queue showing:

- participant and release;
- submitted file and parsed rows;
- automatic field-level score;
- the participant's review note;
- validation or parsing warnings;
- previous submissions and prior decisions.

An administrator may uphold or override the automatic result. Every action
must record the actor, time, old result, new result, and reason. Re-scoring a
release must not silently erase a manual decision.

## Access and certificates

- Anyone with a supported Google, Microsoft, or email-link account may use
  practice releases.
- Practice remains available outside challenge windows.
- Challenge downloads and submissions require enrolment in the applicable
  published round and must obey its opening and closing times.
- A certificate is issued only after the complete Challenge 2 proficiency rule
  is satisfied, including any required manual review.
- The certificate QR code must resolve to a public GHRU Puzzles verification
  page containing a non-guessable certificate code, participant name,
  achievement, issue date, and current validity/revocation state. It must not
  expose private submissions or answer keys.

## Acceptance test

Test the complete participant journey against the uploaded contract 2.1
releases:

1. Register all eight releases from R2.
2. Sign in as a new participant.
3. Download every advertised file, including a zero-byte file.
4. Confirm that deliberately absent files are not rendered as broken links.
5. Download a CSV template and submit all four practice exercises.
6. Confirm conditional `PASS`/`FAIL` validation and scoring.
7. Submit arbitrary outbreak cluster labels and verify partition scoring.
8. Submit blank and `-` for the unavailable typing `wzi` and verify the chosen
   normalisation policy.
9. Open a challenge round, enrol, download, submit, and replace a submission.
10. Close the round and confirm that further downloads/submissions follow the
    configured policy.
11. Request manual review, override a result, and inspect the audit trail.
12. Issue a certificate and verify its QR URL in a signed-out browser.
13. Confirm that no public or participant API can retrieve any object beneath
    a release's `private/` prefix.

The website is ready for Challenge 2 only when this journey succeeds without
database edits, manual object-key construction, or special-case exercise code.
