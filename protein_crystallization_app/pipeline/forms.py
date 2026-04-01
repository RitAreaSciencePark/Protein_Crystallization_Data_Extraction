from django import forms


class SequenceForm(forms.Form):
    sequence = forms.CharField(
        widget=forms.Textarea(attrs={'rows':5, 'cols':50}),
        label="Enter Sequence"
    )
    seq_type_name = forms.CharField(
        max_length=100,
        label="Descriptive Sequence Name"
    )