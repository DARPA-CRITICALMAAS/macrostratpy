from macrostratpy import macrostrat_from_bounds
import argparse


def macrostrat_cli():
    print("Hello from command line")
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--north", type=float, required=True)
    parser.add_argument("-e", "--east", type=float, required=True)
    parser.add_argument("-s", "--south", type=float, required=True)
    parser.add_argument("-w", "--west", type=float, required=True)
    parser.add_argument("-z", "--zoom", type=int,
                        choices=[1, 2, 3, 4, 5, 6, 7, 8, 9], required=True)

    args = parser.parse_args()

    bounds = {
        'n': args.north,
        'e': args.east,
        's': args.south,
        'w': args.west,
    }

    output_path = "/home/efvega/data/macrostrat/hackathon/macrostrat_upper_midwest.json"

    macrostrat_from_bounds(bounds=bounds, output_path=output_path, zoom_level=args.zoom)

    print(args)


if __name__ == "__main__":
    macrostrat_cli()