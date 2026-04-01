from django.shortcuts import render
from django.conf import settings
from .forms import SequenceForm
from .utils import run_pipeline
import os

def run_pipeline_view(request):
    context = {}
    if request.method == "POST":
        form = SequenceForm(request.POST)
        if form.is_valid():
            sequence = form.cleaned_data['sequence']
            seq_type_name = form.cleaned_data['seq_type_name']

            # Run pipeline with MEDIA_ROOT
            result_files = run_pipeline(sequence, seq_type_name, settings.MEDIA_ROOT)

            if result_files is None:
                context['error'] = "No PDB entries found for this sequence."
                context['form'] = form
            else:
                # Prepare relative paths for template
                context['csv_file'] = f"{seq_type_name}/{os.path.basename(result_files['csv'])}"
                context['fasta_file'] = f"{seq_type_name}/{os.path.basename(result_files['fasta'])}"
                context['png_files'] = [f"{seq_type_name}/{os.path.basename(f)}" for f in result_files['pngs']]
                context['pdf_files'] = [f"{seq_type_name}/{os.path.basename(f)}" for f in result_files['pdfs']]
                context['form'] = form
    else:
        form = SequenceForm()
        context['form'] = form

    return render(request, "pipeline/form.html", context)