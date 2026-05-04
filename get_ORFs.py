import argparse
import os
import sys
import requests

from process_ORFs import *

BASE_URL = "https://www.ebi.ac.uk/metagenomics/api/v1"


def fetch_analyses(session, base_url, experiment_type, page_size, page, file_type):
    resp = session.get(
        f"{base_url}/analyses",
        params={
            "page_size": page_size,
            "experiment_type": experiment_type,
            "page": page,
            "file_type": file_type
        },
    )
    resp.raise_for_status()
    data = resp.json().get("data", [])
    return [a["id"] for a in data if "id" in a]


def fetch_downloads(session, base_url, analysis_id):
    resp = session.get(f"{base_url}/analyses/{analysis_id}/downloads")
    resp.raise_for_status()
    return resp.json().get("data", [])

def pick_download(downloads, file_type):

    keywords = {
        'interpro': ['interpro.tsv', 'interpro'],
        'kegg': ['kegg_analysis.tsv', 'kegg_modules'],
        'go': ['_go.tsv', 'go_summary'],
        'cds': ['predicted cds', 'predicted orf', 'predicted_cds', 'predicted_orf']
    }

    targets = keywords.get(file_type, [])

    
    for d in downloads:
        did = (d.get("id") or "").lower()
        desc_obj = d.get("attributes", {}).get("description") or {}
        desc = (desc_obj.get("label", "") + " " + desc_obj.get("description", "")).lower()
        desc = desc.lower()
        if any(t in did or t in desc for t in targets):
            return d
        
    return None

def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--experiment-type",
        choices=["metagenomic", "metatranscriptomic"],
        default="metatranscriptomic",
    )
    p.add_argument(
        "--file-type",
        choices=["go", "kegg", "cds", "interpro"],
        required=True
    )
    p.add_argument("--target-dir", required=True)
    p.add_argument("--page", type=int, default=1)
    p.add_argument("--page-size", type=int, default=500)
    p.add_argument("--max", type=int, default=1000)
    p.add_argument("--base-url", default=BASE_URL)
    args = p.parse_args()

    os.makedirs(args.target_dir, exist_ok=True)

    session = requests.Session()
    processed = 0
    page = args.page

    while processed < args.max:
        analysis_ids = fetch_analyses(
            session=session,
            base_url=args.base_url,
            experiment_type=args.experiment_type,
            page_size=args.page_size,
            page=page,
            file_type=args.file_type
        )
        if not analysis_ids:
            break

        for analysis_id in analysis_ids:
            if processed >= args.max:
                break

            try:
                print(f"fetching downloads for analysis {analysis_id}")

                downloads = fetch_downloads(session, args.base_url, analysis_id)
                chosen = pick_download(downloads, args.file_type)

                if not chosen:
                    print(f"no relevant file in analysis, skipping...")
                    continue

                file_url = chosen.get("links", {}).get("self")
                if not file_url:
                    print(f"no download link for {chosen.get('id')}, skipping…")
                    continue

                out_path = os.path.join(args.target_dir, os.path.basename(chosen["id"]))
                counts_path = f"{os.path.splitext(out_path)[0]}_counts.csv"
                if os.path.exists(counts_path):
                    print("counts already exist, skipping...")

                print(f"downloading {chosen['id']}")
                r = session.get(file_url)
                r.raise_for_status()

                with open(out_path, "wb") as f:
                    f.write(r.content)

                processed += 1

                process(out_path)
                os.remove(out_path)

            except Exception:
                print(f"error with {analysis_id}")
                continue

            
            
        page += 1
        print(processed)
    
    return 0

if __name__ == "__main__":
    raise SystemExit(main())


                


    