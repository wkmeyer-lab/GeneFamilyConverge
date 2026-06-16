"""
convgeno.manifests.checksums — Generate and verify file checksums.

Functions to implement:
    md5sum(path) -> str
        Compute MD5 hash of a file.

    generate_checksum_manifest(file_list, output_path)
        Write a TSV of (filepath, md5, size_bytes, modified_time)
        for all files in the list.

    verify_checksum_manifest(manifest_path) -> list[str]
        Check that all files in a manifest still match their recorded checksums.
        Returns list of mismatched files.

Purpose:
    Reproducibility. "We recorded MD5 checksums of all input files"
    is a one-sentence methods paragraph that reviewers notice.
"""
