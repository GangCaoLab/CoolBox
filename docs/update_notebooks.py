#!/usr/bin/env python
"""
Update all Jupyter notebooks in the docs directory by executing them.
This script runs all notebooks and updates their outputs.
"""

import sys
import argparse
from pathlib import Path
import subprocess
import concurrent.futures
from typing import List, Tuple


def execute_notebook(notebook_path: Path, timeout: int = 300) -> Tuple[Path, bool, str]:
    """
    Execute a single notebook and update its outputs.

    Parameters
    ----------
    notebook_path : Path
        Path to the notebook file
    timeout : int
        Timeout in seconds for notebook execution

    Returns
    -------
    tuple
        (notebook_path, success, message)
    """
    try:
        print(f"Executing: {notebook_path}")

        # Use jupyter nbconvert to execute and update in place
        result = subprocess.run(
            [
                "jupyter", "nbconvert",
                "--to", "notebook",
                "--execute",
                "--inplace",
                "--ExecutePreprocessor.timeout={}".format(timeout),
                "--ExecutePreprocessor.kernel_name=python3",
                str(notebook_path)
            ],
            capture_output=True,
            text=True,
            timeout=timeout + 10  # Add buffer to subprocess timeout
        )

        if result.returncode == 0:
            print(f"✓ Success: {notebook_path}")
            return (notebook_path, True, "Executed successfully")
        else:
            error_msg = result.stderr or result.stdout
            print(f"✗ Failed: {notebook_path}")
            print(f"  Error: {error_msg[:200]}")
            return (notebook_path, False, error_msg)

    except subprocess.TimeoutExpired:
        msg = f"Timeout after {timeout} seconds"
        print(f"✗ Timeout: {notebook_path}")
        return (notebook_path, False, msg)
    except Exception as e:
        msg = str(e)
        print(f"✗ Error: {notebook_path} - {msg}")
        return (notebook_path, False, msg)


def find_notebooks(source_dir: Path, exclude_checkpoints: bool = True) -> List[Path]:
    """
    Find all Jupyter notebooks in the source directory.

    Parameters
    ----------
    source_dir : Path
        Source directory to search
    exclude_checkpoints : bool
        Whether to exclude .ipynb_checkpoints directories

    Returns
    -------
    list
        List of notebook paths
    """
    notebooks = []
    for nb in source_dir.rglob("*.ipynb"):
        if exclude_checkpoints and ".ipynb_checkpoints" in str(nb):
            continue
        notebooks.append(nb)
    return sorted(notebooks)


def main():
    parser = argparse.ArgumentParser(
        description="Execute and update Jupyter notebooks in docs"
    )
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=Path("source"),
        help="Source directory containing notebooks (default: source)"
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=300,
        help="Timeout per notebook in seconds (default: 300)"
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1)"
    )
    parser.add_argument(
        "--filter",
        type=str,
        default="",
        help="Filter notebooks by name pattern (e.g., 'quick_start')"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List notebooks without executing them"
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop on first failure"
    )

    args = parser.parse_args()

    # Find all notebooks
    print(f"Searching for notebooks in {args.source_dir}...")
    notebooks = find_notebooks(args.source_dir)

    # Apply filter
    if args.filter:
        notebooks = [nb for nb in notebooks if args.filter in str(nb)]

    if not notebooks:
        print("No notebooks found.")
        return 0

    print(f"\nFound {len(notebooks)} notebook(s):")
    for nb in notebooks:
        print(f"  - {nb.relative_to(args.source_dir.parent)}")

    if args.dry_run:
        print("\nDry run - no notebooks executed.")
        return 0

    print(f"\nExecuting notebooks (timeout: {args.timeout}s, workers: {args.parallel})...")
    print("=" * 80)

    # Execute notebooks
    results = []
    if args.parallel > 1:
        # Parallel execution
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.parallel) as executor:
            futures = {
                executor.submit(execute_notebook, nb, args.timeout): nb
                for nb in notebooks
            }

            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                results.append(result)

                if args.fail_fast and not result[1]:
                    # Cancel remaining futures
                    for f in futures:
                        f.cancel()
                    break
    else:
        # Sequential execution
        for nb in notebooks:
            result = execute_notebook(nb, args.timeout)
            results.append(result)

            if args.fail_fast and not result[1]:
                break

    # Print summary
    print("\n" + "=" * 80)
    print("Summary:")
    print("=" * 80)

    success_count = sum(1 for _, success, _ in results if success)
    fail_count = len(results) - success_count

    print(f"Total: {len(results)}")
    print(f"Success: {success_count}")
    print(f"Failed: {fail_count}")

    if fail_count > 0:
        print("\nFailed notebooks:")
        for nb, success, msg in results:
            if not success:
                print(f"  ✗ {nb.relative_to(args.source_dir.parent)}")
                print(f"    {msg[:100]}")
        return 1

    print("\n✓ All notebooks executed successfully!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
