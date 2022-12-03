from argparse import ArgumentParser
import logging

from tornado.ioloop import IOLoop
from tornado.log import enable_pretty_logging
from tornado.web import Application

from handlers.version_handler import VersionHandler
from handlers.notfound_handler import NotFoundHandler
from handlers.genome_handler import GenomeHandler
from handlers.query_handler import QueryEngine
from config import log_level

logger = logging.getLogger(__name__)
logger.setLevel(log_level)


def make_app(debug: bool = False) -> Application:
    """Create the Tornado web application.

    Please see https://www.tornadoweb.org/en/stable/guide.html for a crash
    course in how to build web applications with Tornado. Routes are defined
    with regular expressions (the guide contains some simple examples such as
    found at https://tornado-doc-chs.readthedocs.io/en/latest/guide/structure.html#the-application-object).

    """
    logger.info("Starting application")
    return Application(
        [
            ("/version/?", VersionHandler),
            ("/genomehandler/?", GenomeHandler),
            ("/queryengine/(listgenomes)/?", QueryEngine),
            ("/queryengine/(length)/?", QueryEngine),
            ("/queryengine/(retrieveseq)/?", QueryEngine),
            (
                "/queryengine/(searchseq)/?",
                QueryEngine,
            ),
        ],
        debug=debug,
        default_handler_class=NotFoundHandler,
    )


def main() -> None:
    """Run the API service."""
    parser = ArgumentParser()
    parser.add_argument("--debug", "-d", action="store_true", help="enable debug mode")
    parser.add_argument(
        "--port", "-p", type=int, default=8080, help="port to listen on"
    )
    args = parser.parse_args()
    enable_pretty_logging()
    app = make_app(args.debug)
    app.listen(args.port)
    logger.info("Listening on port %d", args.port)
    IOLoop.current().start()


if __name__ == "__main__":
    main()
