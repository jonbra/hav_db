# Database schema normalization proposal

This document describes a proposed normalized schema for the HAV/viral sequence database. It is a proposal only and is NOT applied automatically.

Goals
- Reduce redundancy between sequences and metadata
- Enforce referential integrity with foreign keys
- Make sample-level and sequence-level data clearly separable

High-level entities
- `samples` — one row per biological sample (site, date, patient metadata)
- `submissions` — ingestion events, source file, user, timestamp
- `sequences` — one row per sequence, links to `samples` and `submissions`
- `alignments` — store alignment metadata and path references
- `trees` — store tree builds and newick strings

See `sql/schema_proposal.sql` for a concrete CREATE TABLE proposal and `sql/constraints.sql` for recommended indexes and constraints.
