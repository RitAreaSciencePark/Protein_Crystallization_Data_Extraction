"""
FAIR-Compliant Metadata Generator for Crystal Explorer run outputs
====================================================================

Generates comprehensive metadata for one run's output files according to
FAIR principles:
- Findable: Unique identifiers and rich metadata
- Accessible: Standard formats and clear provenance
- Interoperable: JSON/RDF-compatible structure
- Reusable: Complete documentation and licensing info

Adapted from src/metadata_generator.py (the CLI pipeline's version) to
match Crystal Explorer's own output file names (Output.csv,
Output_compounds.csv, Grouped_conditions.csv, <name>.fasta, ...) instead
of the CLI's (rcsb_hits.csv, merged_results.csv, sequence.fasta, ...).

Author: Nana Njantang Ruth
ORCID: 0000-0002-6003-7521
License: MIT © 2025 RitAreaSciencePark
"""

import os
import json
import hashlib
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any


class FAIRMetadataGenerator:
    """
    Generates FAIR-compliant metadata for one Crystal Explorer run's
    output files.

    Attributes:
        organization_id (str): Unique identifier for RitAreaSciencePark
        package_version (str): Current version of the software
        package_name (str): Name of the software package
    """

    # Organization metadata
    ORGANIZATION_ID = "rit-area-science-park"
    ORGANIZATION_NAME = "RitAreaSciencePark"
    ORGANIZATION_URL = "https://github.com/RitAreaSciencePark"

    # Package metadata
    PACKAGE_NAME = "Protein Crystallization Data Extraction (PCDE) — Crystal Explorer"
    PACKAGE_VERSION = "1.0.0"
    PACKAGE_URL = "https://github.com/RitAreaSciencePark/Protein_Crystallization_Data_Extraction"

    # Author metadata
    AUTHOR_NAME = "Nana Njantang Ruth"
    AUTHOR_ORCID = "0000-0002-6003-7521"

    # License metadata
    LICENSE_TYPE = "MIT"
    LICENSE_URL = "https://opensource.org/licenses/MIT"
    LICENSE_YEAR = 2025

    # File type mappings -- matched against each output file's name with
    # str.endswith(), so more specific suffixes must come before the more
    # general ones they'd otherwise also match (e.g. "Output_compounds.csv"
    # before the plain ".csv"-ish entries).
    FILE_MAPPINGS = {
        "Cryst_cocktail_Table.pdf": "report",
        "Output_compounds.csv": "compounds",
        "Grouped_conditions.csv": "grouped",
        "Output.csv": "hits",
        "PEG.png": "plot_peg",
        "TEMP.png": "plot_temp",
        ".fasta": "fasta",
    }

    def __init__(self, run_identifier: str, output_dir: str, sequence_type: str = None):
        """
        Initialize the metadata generator.

        Args:
            run_identifier (str): The run's folder/protein name -- used to
                name the metadata file and label the dataset, matching
                Crystal Explorer's other per-run output file names.
            output_dir (str): Path to the run's output directory.
            sequence_type (str): Biological sequence type searched
                ('protein', 'dna', 'rna'). Falls back to `run_identifier`
                if not given, so this class still works standalone.
        """
        self.sequence_type_name = run_identifier
        self.sequence_type = sequence_type or run_identifier
        self.output_dir = output_dir
        self.search_timestamp = datetime.utcnow()

    @staticmethod
    def _generate_file_hash(file_path: str, algorithm: str = "sha256") -> str:
        """
        Generate a cryptographic hash of a file for integrity verification.

        Args:
            file_path (str): Path to the file
            algorithm (str): Hash algorithm ('sha256', 'md5')

        Returns:
            str: Hexadecimal hash string
        """
        hash_obj = hashlib.new(algorithm)
        try:
            with open(file_path, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_obj.update(chunk)
            return hash_obj.hexdigest()
        except Exception as e:
            return f"Error computing hash: {str(e)}"

    @staticmethod
    def _get_file_size_kb(file_path: str) -> float:
        """Get file size in kilobytes."""
        try:
            return round(os.path.getsize(file_path) / 1024, 2)
        except Exception:
            return 0.0

    @staticmethod
    def _get_file_mime_type(file_path: str) -> str:
        """Determine MIME type based on file extension."""
        mime_types = {
            ".csv": "text/csv",
            ".pdf": "application/pdf",
            ".png": "image/png",
            ".jpg": "image/jpeg",
            ".jpeg": "image/jpeg",
            ".fasta": "application/x-fasta",
            ".fa": "application/x-fasta",
            ".json": "application/json",
            ".txt": "text/plain",
        }
        ext = Path(file_path).suffix.lower()
        return mime_types.get(ext, "application/octet-stream")

    def _get_file_category(self, file_name: str) -> str:
        """Determine file category from file name."""
        for suffix, category in self.FILE_MAPPINGS.items():
            if file_name.endswith(suffix):
                return category
        return "unknown"

    def _get_file_description(self, file_category: str, file_name: str) -> str:
        """Get a description based on file category and name."""
        descriptions = {
            "fasta": f"Input sequence in FASTA format ({file_name})",
            "hits": f"RCSB sequence search results with X-ray crystallography filters, "
                    f"filtered to entries with experimental crystallization data ({file_name})",
            "compounds": f"Search results with crystallization reagents extracted from "
                         f"pdbx_details ({file_name})",
            "grouped": f"Compounds grouped into unique crystallization conditions, "
                       f"with PDB IDs sharing a condition merged ({file_name})",
            "report": f"Consolidated PDF report with crystallization cocktail details ({file_name})",
            "plot_temp": f"Scatter plot visualization: Temperature analysis ({file_name})",
            "plot_peg": f"Scatter plot visualization: PEG concentration analysis ({file_name})",
            "metadata": f"FAIR-compliant metadata for all output files ({file_name})",
        }
        return descriptions.get(file_category, f"Output file from Crystal Explorer ({file_name})")

    def _generate_file_id(self, file_name: str) -> str:
        """
        Generate a unique file ID using the file name and timestamp.

        Example:
        Grouped_conditions-20260611-154230
        """
        stem = Path(file_name).stem.replace(" ", "_")
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

        return f"{stem}-{timestamp}"

    def generate_file_metadata(self, file_path: str) -> Dict[str, Any]:
        """
        Generate comprehensive metadata for a single output file.

        Args:
            file_path (str): Full path to the file

        Returns:
            Dict: Complete metadata dictionary for the file
        """
        file_name = os.path.basename(file_path)
        file_category = self._get_file_category(file_name)

        # Generate file-specific metadata
        file_metadata = {
            "@context": "https://www.w3.org/ns/dcat",
            "@type": "dcat:Dataset",

            # Unique Identifier (FAIR Principle: Findable)
            "file_id": self._generate_file_id(file_name),
            "file_name": file_name,
            "file_path": os.path.abspath(file_path),
            "relative_path": os.path.relpath(file_path, self.output_dir),

            # File Properties
            "file_properties": {
                "size_kb": self._get_file_size_kb(file_path),
                "file_extension": Path(file_path).suffix,
                "mime_type": self._get_file_mime_type(file_path),
                "encoding": "utf-8" if file_name.endswith(('.csv', '.fasta', '.txt', '.json')) else "binary",
            },

            # Integrity & Validation (FAIR Principle: Accessible)
            "checksum": {
                "algorithm": "SHA-256",
                "value": self._generate_file_hash(file_path, "sha256"),
            },

            # File Classification
            "description": self._get_file_description(file_category, file_name),
            "category": file_category,

            # Input Context
            "input_context": {
                "sequence_type": self.sequence_type,
            },

        }

        return file_metadata

    def generate_dataset_metadata(self, total_files: int = 0, total_size_kb: float = 0) -> Dict[str, Any]:
        """
        Generate metadata for the entire dataset collection.

        Args:
            total_files (int): Total number of files
            total_size_kb (float): Total size in kilobytes

        Returns:
            Dict: Complete dataset metadata
        """
        dataset_metadata = {
            "@context": "https://www.w3.org/ns/dcat",
            "@type": "dcat:Catalog",

            # Dataset Identifier & General Info
            "dataset_id": self.sequence_type_name,
            "dataset_name": f"Crystal Explorer Results - {self.sequence_type_name}",
            "dataset_description": f"Protein crystallization data extraction results for run: "
                                    f"{self.sequence_type_name} (sequence type: {self.sequence_type})",

            # Temporal Information
            "creation_date": self.search_timestamp.isoformat() + "Z",
            "creation_date_human": self.search_timestamp.strftime("%Y-%m-%d %H:%M:%S UTC"),

            # Software Version
            "software": {
                "name": self.PACKAGE_NAME,
                "version": self.PACKAGE_VERSION,
                "url": self.PACKAGE_URL,
            },

            # Organization & Author
            "organization": {
                "id": self.ORGANIZATION_ID,
                "name": self.ORGANIZATION_NAME,
                "url": self.ORGANIZATION_URL,
            },
            "author": {
                "name": self.AUTHOR_NAME,
                "orcid": f"https://orcid.org/{self.AUTHOR_ORCID}",
            },
            "license": {
                "type": self.LICENSE_TYPE,
                "url": self.LICENSE_URL,
                "year": self.LICENSE_YEAR,
            },

            # Data Source Information
            "data_source": {
                "primary_database": "RCSB Protein Data Bank (PDB)",
                "database_url": "https://www.rcsb.org/",
                "search_scope": "Protein structures with X-ray crystallography experimental conditions",
                "api_endpoints": [
                    "https://search.rcsb.org/rcsbsearch/v2/query",
                    "https://data.rcsb.org/rest/v1/core/entry",
                ],
            },

            # Metadata Standard Compliance
            "metadata_standards": {
                "standards": ["FAIR", "Dublin Core", "RDF/JSON-LD"],
                "version": "1.0",
                "compliance_level": "Full compliance with FAIR principles",
            },

            # Access & Distribution
            "access": {
                "url": self.PACKAGE_URL,
                "type": "public",
                "access_level": "open",
                "distribution_format": ["CSV", "FASTA", "PDF", "PNG", "JSON"],
            },

            # Dataset Statistics
            "statistics": {
                "total_files": total_files,
                "total_size_kb": round(total_size_kb, 2),
            },
        }

        return dataset_metadata

    def save_metadata_json(self, file_metadata_list: List[Dict[str, Any]]) -> str:
        """
        Save complete metadata collection to JSON file in the output directory.

        Args:
            file_metadata_list (List[Dict]): List of individual file metadata dictionaries

        Returns:
            str: Path to the saved metadata JSON file
        """
        # Calculate totals
        total_size_kb = round(sum(f.get("file_properties", {}).get("size_kb", 0) for f in file_metadata_list), 2)

        # Generate dataset-level metadata
        dataset_metadata = self.generate_dataset_metadata(
            total_files=len(file_metadata_list),
            total_size_kb=total_size_kb,
        )

        # Compile complete metadata collection
        complete_metadata = {
            "metadata_version": "1.0",
            "metadata_standard": "FAIR Principles Compliant",
            "generation_timestamp": self.search_timestamp.isoformat() + "Z",
            "generation_timestamp_human": self.search_timestamp.strftime("%Y-%m-%d %H:%M:%S UTC"),

            "dataset": dataset_metadata,
            "files": file_metadata_list,

            "summary": {
                "total_files": len(file_metadata_list),
                "total_size_kb": total_size_kb,
                "generation_tool": self.PACKAGE_NAME,
                "generation_tool_version": self.PACKAGE_VERSION,
                "run_identifier": self.sequence_type_name,
                "sequence_type": self.sequence_type,
            },
        }

        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)

        # Save metadata JSON file
        metadata_file_path = os.path.join(self.output_dir, f"{self.sequence_type_name}_metadata.json")

        with open(metadata_file_path, "w", encoding="utf-8") as f:
            json.dump(complete_metadata, f, indent=2, ensure_ascii=False)

        print(f"FAIR metadata saved: {metadata_file_path}")
        return metadata_file_path

    def process_all_output_files(self) -> str:
        """
        Automatically discover and generate metadata for all output files in the output directory.

        Returns:
            str: Path to the generated metadata JSON file
        """
        file_metadata_list = []

        print("Scanning output directory for files...")

        # Walk through the output directory and find all files
        for root, dirs, files in os.walk(self.output_dir):
            for file_name in files:
                # Skip the metadata file itself to avoid circular references
                if file_name.endswith("_metadata.json"):
                    continue

                file_path = os.path.join(root, file_name)

                try:
                    metadata = self.generate_file_metadata(file_path)
                    file_metadata_list.append(metadata)

                    # Show relative path for clarity
                    rel_path = os.path.relpath(file_path, self.output_dir)
                    print(f"  Generated metadata for: {rel_path}")
                except Exception as e:
                    print(f"  Error processing {file_name}: {str(e)}")

        if not file_metadata_list:
            print("  No output files found in the output directory")

        # Save all metadata to a single JSON file
        metadata_json_path = self.save_metadata_json(file_metadata_list)
        return metadata_json_path


def generate_metadata_for_outputs(run_identifier: str, output_dir: str, sequence_type: str = None) -> str:
    """
    Convenience function to generate FAIR-compliant metadata for all output
    files of one Crystal Explorer run.

    Args:
        run_identifier (str): The run's folder/protein name.
        output_dir (str): Path to the run's output directory.
        sequence_type (str): Biological sequence type searched
            ('protein', 'dna', 'rna'). Falls back to `run_identifier` if
            not given.

    Returns:
        str: Path to the generated metadata JSON file

    Example:
        >>> metadata_path = generate_metadata_for_outputs(
        ...     "MyProtein_20260901-120000", "pipeline_outputs/MyProtein_20260901-120000/",
        ...     sequence_type="protein")
        >>> print(f"Metadata saved to: {metadata_path}")
    """
    print("Generating FAIR-compliant metadata for all output files...")
    generator = FAIRMetadataGenerator(run_identifier, output_dir, sequence_type=sequence_type)
    return generator.process_all_output_files()
