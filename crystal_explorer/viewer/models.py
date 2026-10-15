from django.db import models


class SearchRun(models.Model):
    """One record per search submitted from the form -- kept so labs can
    browse past runs (params + outcome) and jump straight back to a
    result page instead of re-running the pipeline."""

    STATUS_RUNNING = "running"
    STATUS_COMPLETED = "completed"
    STATUS_FAILED = "failed"
    STATUS_CHOICES = [
        (STATUS_RUNNING, "Running"),
        (STATUS_COMPLETED, "Completed"),
        (STATUS_FAILED, "Failed"),
    ]

    protein_name = models.CharField(max_length=200)
    folder_name = models.CharField(
        max_length=220, unique=True,
        help_text="Output directory name under PIPELINE_OUTPUT_DIR; also the results/download URL key.",
    )
    sequence = models.TextField(blank=True, help_text="Original sequence used for the search.")
    search_signature = models.CharField(
        max_length=128,
        unique=False,
        db_index=True,
        blank=True,
        help_text="Canonical signature for deduplicating repeated searches.",
    )
    sequence_signature = models.CharField(
        max_length=128,
        unique=False,
        db_index=True,
        blank=True,
        help_text=(
            "Signature of protein name + sequence + sequence type only "
            "(no thresholds/flags) -- a rerun with the same subject but "
            "different parameters reuses and overwrites this row instead "
            "of creating a new history entry."
        ),
    )
    session_key = models.CharField(
        max_length=40, blank=True, default="", db_index=True,
        help_text=(
            "Django anonymous session key of whoever (most recently) ran "
            "this search -- not an account, just a per-browser cookie. "
            "Used to scope the history list to the current visitor and to "
            "gate deletion to the browser that ran the search; the "
            "dedup/reuse logic in index() still checks every row "
            "regardless of this field, so a repeated search is reused "
            "across visitors even though each visitor's history list "
            "shows only their own rows."
        ),
    )
    lab_name = models.CharField(
        max_length=200, blank=True, default="",
        help_text=(
            "Self-reported lab/institution name for this run, typed into "
            "the search form -- optional, never inferred from IP or "
            "session. Surfaced in the run's FAIR metadata JSON."
        ),
    )
    sequence_type = models.CharField(max_length=10)
    identity = models.FloatField()
    evalue = models.FloatField()
    max_hits = models.IntegerField()
    llm_fallback = models.BooleanField(default=False)
    sequence_preview = models.CharField(
        max_length=80, blank=True,
        help_text="First few residues, just enough to recognize the run at a glance.",
    )

    status = models.CharField(max_length=10, choices=STATUS_CHOICES, default=STATUS_RUNNING)
    progress_percent = models.IntegerField(
        default=0,
        help_text="Coarse progress through the pipeline's stages while status=running; 0-100.",
    )
    progress_label = models.CharField(
        max_length=80, blank=True, default="",
        help_text="Short human label for the current stage, e.g. 'Searching RCSB for homologs'.",
    )
    error_message = models.TextField(blank=True)
    row_count = models.IntegerField(null=True, blank=True)
    runtime_seconds = models.FloatField(default=0.0, help_text="Time spent running the pipeline.")

    created_at = models.DateTimeField(auto_now_add=True)
    finished_at = models.DateTimeField(null=True, blank=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self):
        return f"{self.protein_name} ({self.status}) @ {self.created_at:%Y-%m-%d %H:%M}"
