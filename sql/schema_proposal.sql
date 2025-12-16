-- Schema proposal (NOT applied)
-- This file contains CREATE TABLE statements for a normalized schema.

PRAGMA foreign_keys = OFF;

BEGIN TRANSACTION;

CREATE TABLE IF NOT EXISTS submissions (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    source_filename TEXT,
    importer TEXT,
    imported_at TEXT DEFAULT (datetime('now')),
    notes TEXT
);

CREATE TABLE IF NOT EXISTS samples (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT UNIQUE,
    collection_date TEXT,
    country TEXT,
    region TEXT,
    host TEXT,
    sex TEXT,
    age INTEGER,
    other_metadata JSON
);

CREATE TABLE IF NOT EXISTS sequences (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id TEXT UNIQUE,
    sample_id INTEGER REFERENCES samples(id) ON DELETE SET NULL,
    submission_id INTEGER REFERENCES submissions(id) ON DELETE SET NULL,
    gene TEXT,
    length INTEGER,
    sequence TEXT,
    created_at TEXT DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS alignments (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT,
    created_at TEXT DEFAULT (datetime('now')),
    method TEXT,
    parameters TEXT,
    fasta_path TEXT,
    notes TEXT
);

CREATE TABLE IF NOT EXISTS trees (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    alignment_id INTEGER REFERENCES alignments(id) ON DELETE CASCADE,
    created_at TEXT DEFAULT (datetime('now')),
    method TEXT,
    newick TEXT,
    notes TEXT
);

COMMIT;

PRAGMA foreign_keys = ON;
