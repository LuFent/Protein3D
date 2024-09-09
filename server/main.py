from server.flask_server import run_server
from server.tools import import_all_algorithms
from server.algorithms import ALGORITMS

def main():
    directory = "server/algorithms"
    algorithms = ALGORITMS
    run_server(*algorithms)


if __name__ == "__main__":
    main()