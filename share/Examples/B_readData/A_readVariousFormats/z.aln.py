var.warnReadNoFile = False
instring = """CLUSTAL W (1.83) multiple sequence alignment


two             ACCGTATGCTATCGTATTAGCGTATCGTATTAGGCTATT-ATGCGTAATCG---------
four            ACCGTATTCT-TAGCGT-AGCTTAGCGGCATTCGGTACG-ATTCGTATTCGGTTAGCTAG
three           ACCGTATCGTATCGTAT--TCGGTACGTAGTATATTACGTATGCTTATTACGTATTATCG
one             ACCGTACGTACGTATACGTACGTAGGCTAGCTTATCGGC-GCGCTATATCG---------
                ******              *                      *    *

two             ------------
four            CATGCTAGATCG
three           ------------
one             ------------
"""

read(instring)
func.dump()

