"""
URL configuration for config project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/6.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.conf import settings
from django.contrib import admin
from django.contrib.staticfiles.views import serve as serve_static
from django.urls import path, include

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('viewer.urls')),
]

if not settings.DEBUG:
    # Serve static files even with DEBUG=False, since this app has no
    # separate web server (nginx/whitenoise) fronting it locally.
    urlpatterns += [
        path('static/<path:path>', serve_static, {'insecure': True}),
    ]
