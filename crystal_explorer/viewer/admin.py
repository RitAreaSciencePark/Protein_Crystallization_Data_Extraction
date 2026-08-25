from django.contrib import admin

from .models import SearchRun


@admin.register(SearchRun)
class SearchRunAdmin(admin.ModelAdmin):
    list_display = ("protein_name", "status", "row_count", "created_at", "folder_name")
    list_filter = ("status", "sequence_type", "llm_fallback")
    search_fields = ("protein_name", "folder_name")
    readonly_fields = [f.name for f in SearchRun._meta.fields]
