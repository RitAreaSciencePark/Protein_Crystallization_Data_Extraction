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
    error_message = models.TextField(blank=True)
    row_count = models.IntegerField(null=True, blank=True)

    created_at = models.DateTimeField(auto_now_add=True)
    finished_at = models.DateTimeField(null=True, blank=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self):
        return f"{self.protein_name} ({self.status}) @ {self.created_at:%Y-%m-%d %H:%M}"
