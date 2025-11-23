#!/usr/bin/env python3
"""
Compare CPU and GPU function implementations line by line.

1. Find matching functions in CPU (.cpp) and GPU (.cu) files
2. Extract only code inside braces (remove signatures)
3. Remove all comments
4. Compare line by line
"""

import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional


class FunctionExtractor:
    """Extract functions from C++ and CUDA files."""

    @staticmethod
    def find_functions(file_path: str) -> Dict[str, str]:
        """
        Find all functions in a file and return dict of {function_name: body}.
        Body includes only code inside the outermost braces.
        """
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        functions = {}
        lines = content.split('\n')
        i = 0

        while i < len(lines):
            line = lines[i].strip()

            # Skip empty lines and preprocessor directives
            if not line or line.startswith('#'):
                i += 1
                continue

            # Look for function definition with opening brace
            # Pattern: [qualifiers] return_type [class::]function_name(...) [const] {
            pattern = r'(?:__device__|__global__|__host__)?\s*(?:inline\s+)?(?:\w+(?:\s*\*+)?(?:::\w+)?)\s+(\w+(?:::\w+)?)\s*\([^)]*\)\s*(?:const)?\s*\{'

            # Check current line and potentially next few lines for opening brace
            search_text = line
            start_line = i

            # Look ahead up to 5 lines for opening brace
            for j in range(i, min(i + 5, len(lines))):
                if '{' in lines[j]:
                    search_text = ' '.join(lines[i:j+1])
                    break

            match = re.search(pattern, search_text)

            if match:
                func_name = match.group(1)
                # Remove class prefix if present
                if '::' in func_name:
                    func_name = func_name.split('::')[-1]

                # Find where opening brace is
                brace_line = i
                for j in range(i, min(i + 5, len(lines))):
                    if '{' in lines[j]:
                        brace_line = j
                        break

                # Find matching closing brace
                brace_count = 0
                body_start = brace_line
                body_end = None

                for j in range(brace_line, len(lines)):
                    for char in lines[j]:
                        if char == '{':
                            brace_count += 1
                        elif char == '}':
                            brace_count -= 1
                            if brace_count == 0:
                                body_end = j
                                break
                    if body_end is not None:
                        break

                if body_end is not None:
                    # Extract body (everything inside braces)
                    body_lines = []

                    # First line: everything after opening brace
                    first_line = lines[body_start]
                    brace_pos = first_line.find('{')
                    if brace_pos != -1:
                        after_brace = first_line[brace_pos + 1:]
                        if after_brace.strip():
                            body_lines.append(after_brace)

                    # Middle lines: full content
                    for j in range(body_start + 1, body_end):
                        body_lines.append(lines[j])

                    # Last line: everything before closing brace
                    if body_end > body_start:
                        last_line = lines[body_end]
                        brace_pos = last_line.rfind('}')
                        if brace_pos != -1:
                            before_brace = last_line[:brace_pos]
                            if before_brace.strip():
                                body_lines.append(before_brace)

                    body = '\n'.join(body_lines)
                    functions[func_name] = body

                    i = body_end + 1
                    continue

            i += 1

        return functions

    @staticmethod
    def remove_comments(code: str) -> str:
        """Remove all comments from code."""
        # Remove single-line comments
        code = re.sub(r'//.*?$', '', code, flags=re.MULTILINE)

        # Remove multi-line comments
        code = re.sub(r'/\*.*?\*/', '', code, flags=re.DOTALL)

        return code

    @staticmethod
    def normalize_line(line: str) -> str:
        """Normalize a line for comparison."""
        # Remove CUDA qualifiers
        line = re.sub(
            r'__device__|__global__|__host__|__shared__|__constant__', '', line)

        # Normalize whitespace
        line = line.strip()

        # Remove trailing closing braces (for one-liner functions)
        line = re.sub(r'\s*}\s*$', '', line)

        # Normalize math functions
        line = re.sub(r'\bfmin\b', 'min', line)
        line = re.sub(r'\bfmax\b', 'max', line)
        line = re.sub(r'\bfabs\b', 'abs', line)

        return line


def compare_functions(cpu_func: str, gpu_func: str, func_name: str, context_lines: int = 2) -> List[Dict]:
    """
    Compare two function bodies line by line.

    Args:
        cpu_func: CPU function body
        gpu_func: GPU function body
        func_name: Name of the function
        context_lines: Number of lines above/below to check for misalignment
    """
    extractor = FunctionExtractor()

    # Remove comments
    cpu_clean = extractor.remove_comments(cpu_func)
    gpu_clean = extractor.remove_comments(gpu_func)

    # Split into lines
    cpu_lines = [l for l in cpu_clean.split('\n')]
    gpu_lines = [l for l in gpu_clean.split('\n')]

    # Normalize and filter empty lines
    cpu_normalized = [extractor.normalize_line(l) for l in cpu_lines]
    gpu_normalized = [extractor.normalize_line(l) for l in gpu_lines]

    cpu_filtered = [l for l in cpu_normalized if l]
    gpu_filtered = [l for l in gpu_normalized if l]

    # Compare line by line with misalignment detection
    differences = []
    max_len = max(len(cpu_filtered), len(gpu_filtered))

    for i in range(max_len):
        cpu_line = cpu_filtered[i] if i < len(cpu_filtered) else None
        gpu_line = gpu_filtered[i] if i < len(gpu_filtered) else None

        if cpu_line != gpu_line:
            # Check if this is a misalignment by looking at neighboring lines
            is_misalignment = False

            if cpu_line and gpu_line:
                # Check if CPU line appears in nearby GPU lines
                for offset in range(-context_lines, context_lines + 1):
                    if offset == 0:
                        continue
                    gpu_idx = i + offset
                    if 0 <= gpu_idx < len(gpu_filtered):
                        if cpu_line == gpu_filtered[gpu_idx]:
                            is_misalignment = True
                            break

                # Check if GPU line appears in nearby CPU lines
                if not is_misalignment:
                    for offset in range(-context_lines, context_lines + 1):
                        if offset == 0:
                            continue
                        cpu_idx = i + offset
                        if 0 <= cpu_idx < len(cpu_filtered):
                            if gpu_line == cpu_filtered[cpu_idx]:
                                is_misalignment = True
                                break

            # Only record if it's not just a misalignment
            if not is_misalignment:
                differences.append({
                    'line_num': i + 1,
                    'cpu': cpu_line if cpu_line else '(missing)',
                    'gpu': gpu_line if gpu_line else '(missing)'
                })

    return differences


def find_file_pairs(base_dir: str) -> List[Tuple[str, str]]:
    """Find matching CPU and GPU file pairs."""
    pairs = []

    cpu_dir = Path(base_dir) / "cpu-shower"
    gpu_dir = Path(base_dir) / "gpu-shower"

    if not cpu_dir.exists() or not gpu_dir.exists():
        print(f"Error: cpu-shower or gpu-shower not found in {base_dir}")
        return pairs

    # Find all .cpp files
    for cpu_file in cpu_dir.rglob("*.cpp"):
        rel_path = cpu_file.relative_to(cpu_dir)
        gpu_file = gpu_dir / str(rel_path).replace('.cpp', '.cu')

        if gpu_file.exists():
            pairs.append((str(cpu_file), str(gpu_file)))

    return pairs


def main():
    """Main comparison routine."""
    base_dir = Path(__file__).parent

    print("=" * 80)
    print("CPU vs GPU Function Comparison")
    print("=" * 80)
    print()

    # Find file pairs
    pairs = find_file_pairs(base_dir)
    print(f"Found {len(pairs)} file pairs")
    print()

    extractor = FunctionExtractor()
    total_functions = 0
    total_matches = 0
    total_diffs = 0

    all_results = []

    # Compare each file pair
    for cpu_file, gpu_file in pairs:
        cpu_name = Path(cpu_file).name
        gpu_name = Path(gpu_file).name

        # Extract functions
        cpu_functions = extractor.find_functions(cpu_file)
        gpu_functions = extractor.find_functions(gpu_file)

        # Find common functions
        common = set(cpu_functions.keys()) & set(gpu_functions.keys())

        if not common:
            continue

        print(f"\n{cpu_name} <-> {gpu_name}")
        print("-" * 80)

        for func_name in sorted(common):
            total_functions += 1

            diffs = compare_functions(
                cpu_functions[func_name],
                gpu_functions[func_name],
                func_name
            )

            if diffs:
                total_diffs += 1
                print(f"\n  ✗ {func_name}: {len(diffs)} differences")

                # Show first few differences
                for diff in diffs[:5]:
                    print(f"    Line {diff['line_num']}:")
                    print(f"      CPU: {diff['cpu']}")
                    print(f"      GPU: {diff['gpu']}")

                if len(diffs) > 5:
                    print(f"    ... and {len(diffs) - 5} more differences")

                all_results.append({
                    'file': cpu_name,
                    'function': func_name,
                    'differences': diffs
                })
            else:
                total_matches += 1
                print(f"  ✓ {func_name}: match")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total functions compared: {total_functions}")
    print(f"Matches: {total_matches}")
    print(f"Differences: {total_diffs}")
    print()

    # Write detailed report
    if all_results:
        report_file = base_dir / "function_comparison.dat"
        with open(report_file, 'w') as f:
            f.write("CPU vs GPU Function Comparison - Detailed Report\n")
            f.write("=" * 80 + "\n\n")

            for result in all_results:
                f.write(f"\nFile: {result['file']}\n")
                f.write(f"Function: {result['function']}\n")
                f.write(f"Differences: {len(result['differences'])}\n")
                f.write("-" * 80 + "\n")

                for diff in result['differences']:
                    f.write(f"\nLine {diff['line_num']}:\n")
                    f.write(f"  CPU: {diff['cpu']}\n")
                    f.write(f"  GPU: {diff['gpu']}\n")

        print(f"Detailed report written to: {report_file}")

    return 0 if total_diffs == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
