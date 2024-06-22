from server.flask_server import run_server
from server.tools import import_all_algorithms


def main():
    directory = "server/algorithms"
    algorithms = import_all_algorithms(directory)
    run_server(*algorithms)


if __name__ == "__main__":
    main()