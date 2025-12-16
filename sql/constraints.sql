-- Recommended constraints, indexes and integrity checks

-- Indexes for faster joins/filters
CREATE INDEX IF NOT EXISTS idx_sequences_sample_id ON sequences(sample_id);
CREATE INDEX IF NOT EXISTS idx_sequences_submission_id ON sequences(submission_id);
CREATE INDEX IF NOT EXISTS idx_samples_sample_id ON samples(sample_id);

-- Example CHECK constraints (SQLite limited support)
-- Ensure length is non-negative
-- Note: CHECKs may be refused by older SQLite versions; review before applying
ALTER TABLE sequences ADD COLUMN verified_length INTEGER DEFAULT 0;

-- Foreign key enforcement is enabled at connection time: PRAGMA foreign_keys = ON;
