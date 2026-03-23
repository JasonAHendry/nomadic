import socket


def next_free_port(port: int, max_port: int = 65535):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    while port <= max_port:
        try:
            sock.bind(("", port))
            sock.close()
            return port
        except OSError:
            port += 1
    raise IOError(f"Could not find a port from {port} to {max_port}")
