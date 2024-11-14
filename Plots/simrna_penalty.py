import matplotlib.pyplot as plt
import matplotlib
import signal

signal.signal(signal.SIGINT, signal.SIG_DFL)


def main():
    matplotlib.rc("font", family="sans-serif")
    matplotlib.rc("font", serif="Helvetica Neue")
    matplotlib.rcParams.update({"font.size": 14})
    plt.rcParams["text.usetex"] = True
    # plt.rcParams.update({"figure.dpi": 300})

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot([0, 3.5, 3.5, 9.5, 9.5, 25, 30], [3.5, 0, -1, -1, 0, 25, 25], lw=3, c="red")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(r"distance [\AA]")
    ax.set_ylabel(r"score")
    ax.set_yticks([0, 3.5])
    ax.set_xticks([0, 3.5, 9.5, 25])
    ax.axhline(y=0, color="k", lw=3)
    ax.axvline(x=0, color="k", lw=3)
    ax.grid()
    return fig, ax


if __name__ == "__main__":
    fig, ax = main()
    fig.tight_layout()
    plt.savefig("penalty.png", dpi=300)
