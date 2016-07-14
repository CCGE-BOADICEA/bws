""" Command line utility. """
import os
from django.core.management.base import BaseCommand, CommandError
from boadicea.pedigree import PedigreeFile


class Command(BaseCommand):
    help = 'Convert old style pedigree format to BOADICEA import pedigree v4'

    def add_arguments(self, parser):
        parser.add_argument('file', type=str)

    def handle(self, *args, **options):
        pedigree_file = options['file']
        if not os.path.isfile(pedigree_file):
            raise CommandError('File "%s" does not exist' % pedigree_file)
        if not os.access(pedigree_file, os.R_OK):
            raise CommandError('File "%s" can not be read' % pedigree_file)

        PedigreeFile.convert_boadice_pedigree_v4(pedigree_file)
