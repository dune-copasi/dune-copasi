#!/usr/bin/env python3

from jschon import create_catalog, JSON, JSONSchema
import pprint
import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Simple CLI for JSCHON", exit_on_error=False)
    parser.add_argument('-q', '--quiet', action='store_true')
    parser.add_argument('-c', '--catalog', default='2020-12')
    parser.add_argument('-s', '--schema')
    parser.add_argument('-j', '--json')
    args = parser.parse_args()

    create_catalog(args.catalog)

    print(args.json)
    json = JSON.loadf(os.path.realpath(args.json))
    schema = JSONSchema.loadf(os.path.realpath(args.schema))
    result = schema.evaluate(json)
    if not args.quiet:
      pprint.pp(result.output('basic'))
    if result.valid:
      return sys.exit(0)
    else:
      return sys.exit(1)

if __name__ == '__main__':
    main()