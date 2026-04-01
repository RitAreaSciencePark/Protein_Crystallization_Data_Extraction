from django.urls import path
from . import views

urlpatterns = [
    path('', views.run_pipeline_view, name='run_pipeline'),
]