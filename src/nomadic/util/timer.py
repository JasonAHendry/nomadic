import time


class Timer:
    """
    Simple timer class for measuring execution time of code blocks
    """

    def __init__(self, name: str):
        self.start_time = None
        self.timings = {}
        self.name = name

    def start(self):
        self.start_time = time.time()

    def time(self, name: str):
        if self.start_time is None:
            raise RuntimeError("time can only be called after start")
        end_time = time.time()
        elapsed_time = end_time - self.start_time
        self.timings[name] = elapsed_time
        self.start_time = end_time

    def report(self):
        print(f"Execution time report {self.name}:")
        for name, timing in self.timings.items():
            print(f"  {name}: {timing:.2f} seconds")
