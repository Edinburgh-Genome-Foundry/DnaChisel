from dnachisel.tailor.SequenceDesigner import SequenceDesigner

from dnachisel.builtin_specifications import EnforceTranslation
from dnachisel.tailor.Specification import CAI,GCContent
from dnachisel.tailor.DnaDesignProblem import DnaDesignProblem

class myDesigner(SequenceDesigner):

    def define_problem(self, sequence,solution_id):
        problem = DnaDesignProblem(
            sequence =sequence,
            solution_id = solution_id,
#             constraints = [
#                 EnforceTranslation(),
#             ],
            design_space = self.design_space
        )
        
        return problem
        


from dnachisel.tailor.Design import DesignSpace,FullFactorial,Optimization
from dnachisel.tailor.tools import Range, RangeSet
from dnachisel.builtin_specifications.EnforceTranslation import EnforceTranslation
from dnachisel.tailor.Specification import CAI, GCContent

# design =   FullFactorial(
design =      Optimization(
        ['cai','gc'],
        [GCContent(),GCContent()],
        [RangeSet([(1, Range(0,0.06)),(3, Range(0.06,0.07)),(2, Range(0.07,1))]),   
         RangeSet([('a', Range(0,0.3)),('b', Range(0.3,0.6)),('c', Range(0.6,1))]) ],
        ['REAL','REAL'],
        target = '2.a'
)

designer = myDesigner(
    "test", 
    'cggcttaactcgagagctgacctgcttctacctagccttacagtggtaacgacccaatctgcgtagcgcaacgca', 
    design, 
    "./test", 
    createDB=True
)

designer.run()