from django.db import models

class Submission(models.Model):
    sequence_name = models.CharField(max_length=100)
    sequence = models.TextField()
    output_file = models.FileField(upload_to='outputs/')
    created_at = models.DateTimeField(auto_now_add=True)