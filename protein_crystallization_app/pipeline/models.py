
from django.db import models


class Submission(models.Model):
    sequence_name = models.CharField(max_length=100)
    sequence = models.TextField()

    # core outputs
    fasta_file = models.FileField(upload_to='outputs/fasta/')
    rcsb_csv = models.FileField(upload_to='outputs/rcsb/')
    merged_csv = models.FileField(upload_to='outputs/merged/')

    # optional visual outputs
    plot_files = models.JSONField(default=list, blank=True)  # list of PNG paths
    pdf_files = models.JSONField(default=list, blank=True)   # list of PDF paths

    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.sequence_name