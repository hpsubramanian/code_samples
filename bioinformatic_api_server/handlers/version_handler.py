from handlers.base_handler import BaseView
import tornado


class VersionHandler(BaseView):
    """Handle version requests."""

    SUPPORTED_METHODS = (
        "GET",
        "POST",
    )

    def get(self) -> None:
        """Write version info."""
        self.write({"api": self.__version__, "tornado": tornado.version})

    def post(self) -> None:
        """Why would you POST to a version endpoint?"""
        self.set_status(405, reason="you can't change the version via POST")
