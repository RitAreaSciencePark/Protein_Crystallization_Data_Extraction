from django.urls import path

from . import views

app_name = "viewer"

urlpatterns = [
    path("", views.home, name="home"),
    path("explorer/", views.index, name="index"),
    path("organization/", views.organization, name="organization"),
    path("history/", views.history, name="history"),
    path("history/<int:run_id>/delete/", views.delete_run, name="delete_run"),
    path("protein/<str:protein_name>/", views.results, name="results"),
    path("protein/<str:protein_name>/download/<str:kind>/", views.download, name="download"),
]
