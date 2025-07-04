import logging


class LoggingFascade:
    """
    Interface with pythons logging module

    """

    # Formatting defaults
    DATE_FORMAT = "%Y-%m-%d %H:%M"
    STREAM_FORMAT = "%(message)s"
    FILE_FORMAT = "[%(asctime)s][%(levelname)s] %(message)s"

    def __init__(
        self, logger_name: str = "Default", verbose: bool = False, log_path: str = None
    ):
        """
        Instantiate the the default logger, create the console and optionally
        the file handler

        """

        self.logger = logging.getLogger(logger_name)
        self.logger.setLevel(logging.DEBUG if verbose else logging.INFO)

        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self._add_console_handler()

        self.log_path = log_path
        if self.log_path is not None:
            self._add_file_handler(log_path)

    def _add_console_handler(self):
        """
        Add a console handler

        """

        console_handler = logging.StreamHandler()
        console_formatter = logging.Formatter(self.STREAM_FORMAT)
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)

    def _add_file_handler(self, log_path: str):
        """
        Add a file handler

        """

        file_handler = logging.FileHandler(log_path)
        file_formatter = logging.Formatter(self.FILE_FORMAT, self.DATE_FORMAT)
        file_handler.setFormatter(file_formatter)
        self.logger.addHandler(file_handler)

    def info(self, message: str):
        self.logger.info(message)

    def debug(self, message: str):
        self.logger.debug(message)

    def warning(self, message: str):
        self.logger.warning(message)

    def error(self, message: str):
        self.logger.error(message)

    def critical(self, message: str):
        self.logger.critical(message)
