from server.server import run_server
from server.tools import import_all_algorithms


def main():
    directory = "server/algorithms"  # Replace with the actual path to your directory
    algorithms = import_all_algorithms(directory)
    print(algorithms)
    run_server(*algorithms)
