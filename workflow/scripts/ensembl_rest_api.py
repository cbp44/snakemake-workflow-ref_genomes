#!/usr/bin/env python3
# A Python interface for accessing the Ensembl REST API
# It is used to get metadata about reference genomes

import argparse
import sys
import json
import time
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError


class EnsemblRestClient:
    def __init__(self, server="https://rest.ensembl.org", reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(
        self, endpoint: str, headers: dict = None, params: dict = None
    ):
        """Perform a REST API call to the given endpoint.

        Args:
            endpoint (str): Server endpoint to perform call
            headers (dict, optional): HTTP headers to send when making call. Defaults to None.
            params (dict, optional): HTTP parameters to add to endpoint for call. Defaults to None.

        Returns:
            [object]: Python object containing JSON document of response or None on error
        """
        if headers is None:
            headers = {}

        if "Content-Type" not in headers:
            headers["Content-Type"] = "application/json"

        if params:
            endpoint += "?" + urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = Request(self.server + endpoint, headers=headers)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if "Retry-After" in e.headers:
                    retry = e.headers["Retry-After"]
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, headers, params)
            else:
                sys.stderr.write(
                    "Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n".format(
                        endpoint, e
                    )
                )

        return data

    def get_genome_info(self, species: str):
        """Get Ensembl genome info for a given species.

        Args:
            species (str): Species to get data for e.g. homo_sapiens

        Returns:
            object: JSON document of genome info for species
        """
        return self.perform_rest_action("/info/genomes/{0}".format(species))

    def get_assembly_info(self, species: str):
        return self.perform_rest_action("/info/assembly/{0}".format(species))

    def get_variation_consequence_types(self, species: str, incl_rank: int = 1):
        return self.perform_rest_action(
            "/info/variation/consequence_types/?rank={0}".format(incl_rank)
        )


def variation_consequence_func(args):
    client = EnsemblRestClient()

    info = client.get_variation_consequence_types(args.species)
    info_pp = json.dumps(info, indent=4, sort_keys=True)
    args.output.write(info_pp)


def genome_info_func(args):
    client = EnsemblRestClient()

    info_pp = json.dumps(client.get_genome_info(args.species), indent=4, sort_keys=True)
    args.output.write(info_pp)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Ensembl REST API interface")
    parser.add_argument("--species", type=str, default="homo_sapiens")

    subparsers = parser.add_subparsers()

    genome_info = subparsers.add_parser("genome_info")
    genome_info.add_argument(
        "--output", default=sys.stdout, type=argparse.FileType("w", encoding="UTF-8")
    )
    genome_info.set_defaults(func=genome_info_func)

    variation_consequence = subparsers.add_parser("variation_consequences")
    variation_consequence.add_argument(
        "--output", default=sys.stdout, type=argparse.FileType("w", encoding="UTF-8")
    )
    variation_consequence.set_defaults(func=variation_consequence_func)

    args_to_parse = sys.argv[1:]
    is_snakemake = False
    try:
        if snakemake:
            is_snakemake = True
            args_to_parse = snakemake.params.extra.split(" ")
    except:
        pass
    args = parser.parse_args(args_to_parse)
    args.func(args)
