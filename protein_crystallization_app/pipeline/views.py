# views.py
from django.conf import settings
from django.shortcuts import render
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

            result_files = run_pipeline(sequence, seq_type_name, settings.MEDIA_ROOT)

            if result_files is None:
                context['error'] = "No PDB entries found for this sequence."
            else:
                csv = result_files.get('merged_csv')
                fasta = result_files.get('fasta')

                context['csv_file'] = (
                    f"{seq_type_name}/{os.path.basename(csv)}"
                    if isinstance(csv, (str, bytes, os.PathLike)) else None
                )

                context['fasta_file'] = (
                    f"{seq_type_name}/{os.path.basename(fasta)}"
                    if isinstance(fasta, (str, bytes, os.PathLike)) else None
                )

                # Ensure pngs is a flat list of valid paths
                pngs = result_files.get('pngs', [])
                context['png_files'] = [
                    f"{seq_type_name}/{os.path.basename(f)}"
                    for f in pngs
                    if isinstance(f, (str, bytes, os.PathLike))
                ]

                # Flatten pdfs in case it's nested
                raw_pdfs = result_files.get('pdfs', [])
                pdfs = []
                for item in raw_pdfs:
                    if isinstance(item, list):
                        pdfs.extend(item)
                    else:
                        pdfs.append(item)

                context['pdf_files'] = [
                    f"{seq_type_name}/{os.path.basename(f)}"
                    for f in pdfs
                    if isinstance(f, (str, bytes, os.PathLike))
                ]

            context['form'] = form
    else:
        context['form'] = SequenceForm()

    return render(request, "pipeline/form.html", context)