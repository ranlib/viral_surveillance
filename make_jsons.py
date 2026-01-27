#!/usr/bin/env python
"""
Purpose:
make json file for a list of fastq files

Input: 
1) ABSOLUTE directory path to fastq files 
2) json file template

Output: 
1) for each sample json file with sample name

Note: sample name created from fastq file name with the assumption
that fastq file name has the structure samplename_XXXXX.fastq.gz

"""
import sys
import os
import sys
import json
import re
import pathlib
import argparse

parser = argparse.ArgumentParser(prog="make_json", description="make input json files for varpipe for fastq files in input directory")
parser.add_argument("--dir", "-d", required=True, type=str, help="Input directory with fastq files")
parser.add_argument("--json_template", "-t", required=True, type=argparse.FileType("r"), help="Output json file name")
parser.add_argument("--json_output_dir", "-o", required=True, type=str, help="Output directory for json files")
parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
args = parser.parse_args()

files = []
for filepath in pathlib.Path(args.dir).glob("**/*.fastq.gz"):
    files.append(str(filepath.absolute()))

if args.verbose:
    print(files)

with open(args.json_template.name, "r", encoding="ascii") as stream:
    try:
        json_template = json.load(stream)
    except json.JSONDecodeError as exc:
        print(exc)

R1s = [r for r in files if "_R1" in r]
if len(R1s) == 0:
    print(f"<E> make_jsons: length of fastq list = {len(R1s)}")
    sys.exit(1)

records = []
for r1 in R1s:
    r2 = r1.replace("_R1", "_R2")
    filename_without_ext = pathlib.Path(r1).stem
    sample_id = filename_without_ext.split('_', 1)[0]
    records.append({"sample_id": sample_id, "fastq_R1": r1, "fastq_R2": r2})
    
json_template["ViralSurveillanceCohort.samples"] = records

output_json = pathlib.Path(args.json_output_dir, "test.json")
with open(output_json, "w", encoding="ascii") as f:
    json.dump(json_template, f, indent=4)
