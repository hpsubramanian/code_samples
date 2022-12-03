from handlers.base_handler import BaseView
import tornado


class NotFoundHandler(BaseView):
    """Default handler for unmapped requests"""

    def prepare(self):
        raise tornado.web.HTTPError(status_code=404, reason="Invalid resource path.")
