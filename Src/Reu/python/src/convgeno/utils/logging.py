"""
convgeno.utils.logging — Configure pipeline logging.

Functions to implement:
    setup_logger(name, log_dir=None, level="INFO") -> logging.Logger
        Create a logger that writes to both console and (optionally) a file.
        Console output is human-readable; file output includes timestamps.

    log_step_start(logger, step_name)
        Log the start of a pipeline step with a clear separator.

    log_step_end(logger, step_name, success=True)
        Log the completion of a pipeline step.
"""
