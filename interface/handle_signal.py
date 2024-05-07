import signal
import sys
import time


def signal_handler(sig, frame):
    print()
    print("SIGNAL RECEIVED: shutting down...")
    print()
    sys.exit(0)

if __name__ == "__main__":

    # Register signal handler function
    signal.signal(signal.SIGINT, signal_handler)

    # Run loop intil ctrl-c is pressed
    while True:
        time.sleep(.5)
        print("running & waiting for signal!")