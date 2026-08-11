from django.urls import path

from . import views

app_name = "viewer"

urlpatterns = [
    path("", views.index, name="index"),
    path("protein/<str:protein_name>/", views.results, name="results"),
    path("protein/<str:protein_name>/download/<str:kind>/", views.download, name="download"),
]
