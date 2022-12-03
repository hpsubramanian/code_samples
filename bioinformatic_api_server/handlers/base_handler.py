from tornado.web import RequestHandler

import json
from config import __version__


class BaseView(RequestHandler):
    """
    Base Handler for this application.
    """

    # https://opensource.com/article/18/6/tornado-framework

    def initialize(self):
        """To set version number, Tornado preferred way"""
        self.__version__ = __version__

    def prepare(self):
        """Convert incoming data to utf8"""
        self.form_data = {
            key: [val.decode("utf8") for val in val_list]
            for key, val_list in self.request.arguments.items()
        }

    def set_default_headers(self):
        """Set the default response header to be JSON."""
        self.set_header("Content-Type", "application/json;")

    def send_response(self, data, status=200):
        """Construct and send a JSON response with appropriate status code."""
        self.set_status(status)
        response_packet = {
            "api_version": self.__version__,
            "data": data,
            "status": status,
        }
        self.write(json.dumps(response_packet))
