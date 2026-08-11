from django import forms


class SequenceSearchForm(forms.Form):
    protein_name = forms.CharField(
        label="Protein name",
        max_length=200,
        widget=forms.TextInput(attrs={"placeholder": "e.g. Lysozyme"}),
    )
    sequence = forms.CharField(
        label="Sequence (raw, or pasted FASTA)",
        widget=forms.Textarea(attrs={
            "rows": 8,
            "placeholder": ">my_protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEV...",
        }),
    )
    sequence_type = forms.ChoiceField(
        label="Sequence type",
        choices=[("protein", "Protein"), ("dna", "DNA"), ("rna", "RNA")],
        initial="protein",
    )
    identity = forms.FloatField(
        label="Minimum identity", initial=0.5, min_value=0.0, max_value=1.0,
        help_text="Fraction, 0-1",
    )
    evalue = forms.FloatField(
        label="Max E-value", initial=1.0, min_value=0.0,
    )
    max_hits = forms.IntegerField(
        label="Max PDB entries", initial=25, min_value=1, max_value=200,
        help_text="Kept modest by default since this runs synchronously in your browser request.",
    )
    llm_fallback = forms.BooleanField(
        label="Use Claude to resolve compounds the dictionary misses",
        required=False,
        help_text="Requires ANTHROPIC_API_KEY set on the server.",
    )

    def clean_sequence(self):
        """Strip FASTA header lines ('>' prefixed) and whitespace, same
        parsing behaviour as the CLI's prompt_for_sequence()/read_fasta()."""
        raw = self.cleaned_data["sequence"]
        lines = [ln.strip() for ln in raw.splitlines()]
        seq_lines = [ln for ln in lines if ln and not ln.startswith(">")]
        sequence = "".join(seq_lines).strip().upper()
        if not sequence:
            raise forms.ValidationError("No sequence found (paste a raw sequence or a FASTA record).")
        return sequence

    def clean_protein_name(self):
        name = self.cleaned_data["protein_name"].strip()
        if not name:
            raise forms.ValidationError("Protein name is required.")
        return name
