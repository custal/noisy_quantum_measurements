(* ::Package:: *)

(*The functions in this file are used to calculate results for arbitrary quantum measurement protocols. In this code protocol could either refer to the assignment of outcomes as a list of numbers or to the complete results for a given assignment*)


(*Trine POVM*)
phi1Trine={{1},{0}};
phi2Trine={{1/2},{Sqrt[3]/2}};
phi3Trine={{1/2},{-Sqrt[3]/2}};

(*Defining the matrices M0, M1,and M2*)
povmsTrine={2/3*Matrixify[phi1Trine],
2/3*Matrixify[phi2Trine],
2/3*Matrixify[phi3Trine]};

outcomesTrine = PovmPermutations[povmsTrine];

(*Tetrahedron POVM*)
phi1Tetra={1,0};
phi2Tetra={Rationalize[1/Sqrt[3]],Rationalize[Sqrt[2/3]]};
phi3Tetra={Rationalize[1/Sqrt[3]],E^(I*2*Pi/3)*Rationalize[Sqrt[2/3]]};
phi4Tetra={Rationalize[1/Sqrt[3]],E^(-I*2*Pi/3)*Rationalize[Sqrt[2/3]]};

povmsTetra = { 1/2*Matrixify[phi1Tetra],
		1/2*Matrixify[phi2Tetra],
		1/2*Matrixify[phi3Tetra],
		1/2*Matrixify[phi4Tetra]};

outcomesTetra = PovmPermutations[povmsTetra];

(*Pre calculated results for optimal protocols*)
GetProtocolsTrine := Import[FileNameJoin[{ParentDirectory[NotebookDirectory[]], "data", "trine03082023.mx"}]];
GetProtocolsTetra := Import[FileNameJoin[{ParentDirectory[NotebookDirectory[]], "data", "tetra03082023.mx"}]];


MeasureAbsolute[x_, y_] := ((1 - x)^2 + y^2)/(1 - x + y)
MeasureAbsolute::usage = 
"
x (Number): probability of measuring 0 on input state |0>
y (Number): probability of measuring 0 on input state |1>

Returns (Number): the figure of merit for an x and y value

This is the figure of merit based of the absolute value.
";

OptimalBasisValues[eigVals_] := 
 Module[{x, y, \[Lambda]min = Min[eigVals], \[Lambda]max = Max[eigVals]},
	  {x, y} = Piecewise[{
	     {{\[Lambda]max, (Sqrt[2] - 1)*(1 - \[Lambda]max)}, \[Lambda]max <=
	        1 - (Sqrt[2] + 1)*\[Lambda]min}, 
	     {{\[Lambda]max, \[Lambda]min}, 
	      1 - (Sqrt[2] + 1)*\[Lambda]min < \[Lambda]max <= 
	       1 - (Sqrt[2] - 1)*\[Lambda]min},
	     {{1 - (Sqrt[2] - 1)*\[Lambda]min, \[Lambda]min}, 
	      1 - (Sqrt[2] - 1)*\[Lambda]min < \[Lambda]max}
	     }];
	  {x, y}
	  ]
OptimalBasisValues::usage =
"
eigVals (List[Number]): list of eigenvalues

Returns (List[Number]): The optimal x and y values based on the minimum and maximum eigenvalues

This is for the figure of merit based of the absolute value.
";
  
PovmPermutations[povms_] := 
 Module[{outcomes = <||>}, 
	  Do[outcomes[{i, j}] = KroneckerProduct[povms[[i]], povms[[j]]];
	   , {i, Length[povms]}, {j, Length[povms]}]; 
	  outcomes]
PovmPermutations::usage =
"
povms (List[Matrix]): list of POVM elements 

Returns (Association[List[Number] -> Matrix]): The KroneckerProduct of all pairwise combinations of the input povms

TODO: Generalise this to an arbitrary number of systems, including 1.
";

Matrixify[v_] := Outer[Times, v, ConjugateTranspose[v]]
Matrixify::usage =
"
v (List[Number]): A one dimensional list representing a ket. Should be of the form {a, b, ...} 

Returns (Matrix): The outer product of the input list v with itself

This puts a pure state into density matrix form.
";

GenerateProtocols[outcomes_] := 
 Module[{permutations, protocols},
  
	  permutations = Keys[outcomes];
	  protocols = <|{} -> Null,
	    permutations -> Null|>;
	  
	  Do[
	   protocols[selection] = Null;
	   protocols[Complement[permutations, selection]] = Null;
	   , {selection, Subsets[permutations, Floor[Length[outcomes]/2]]}];
	  protocols]
GenerateProtocols::usage =
"
outcomes (Association[array]): Assocation of all possible outcomes and their corresponding POVM

Returns (Association[List[number] -> Null]): An association with keys corresponding to all possible combinations of outcomes and Null values

For a set of possible outcomes, this function will generate all the possible groupings of these outcomes and place them as keys in an association with Null entries which will later be populated with numerical results for each assignment.
";

  
SumPovms[protocol_, povmOutcomes_] := Module[{sum},
(*Use the first povm matrix to find the dimension assuming all povms have the same dimension*)
	  sum = ConstantArray[0, Dimensions[Values[povmOutcomes][[1]]]];
	  
	  Do[sum += povmOutcomes[outcome]
	   , {outcome, protocol}];
	  sum]
SumPovms::usage =
"
protocol (List[List[number]]): A list containing the outcomes for an assignment
povmOutcomes (Association[List[Number] -> Matrix]): An association where the keys are the outcomes and the values are the corresponding POVMs

Returns (Matrix): The POVM matrix of the input assignment (protocol)

Sum together the POVM matrices of the input outcomes assignment.
";

CalculateMeasures[protocol_, povmOutcomes_] := 
 Module[{Q, eigenValues, eigenVectors, x, y, measureValue, result},
	  
	  Q = SumPovms[protocol, povmOutcomes];
	  {eigenValues, eigenVectors} = Chop[N[Eigensystem[Q]]];
	  {x, y} = OptimalBasisValues[eigenValues];
	  measureValue = MeasureAbsolute[x, y];
	 
	  result = <|"protocol" -> protocol, "eigenValues" -> eigenValues, 
	    "eigenVectors" -> eigenVectors, "measureValue" -> measureValue, 
	    "measureFunction" -> "MeasureAbsolute", "xStarYStar" -> {x, y}|>;
	  result
	  ]
CalculateMeasures::usage =
"
protocol (List[List[number]]): A list containing the outcomes for an assignment
povmOutcomes (Association[List[Number] -> Matrix]): An association where the keys are the outcomes and the values are the corresponding POVMs

Returns (Association[List[number] -> Association[...]]): An association with keys corresponding to all possible combinations of outcomes and calculated results in the values

Given a protocol, calculate the eigenvalues, eigenvectors, x*, y* and figure of merit ffrom the sum of the protocol POVMs. measureFunction is hard coded as currently there is only one.
";
  
CalculateAllProtocols[outcomes_] := 
 Module[{protocols, protocol, result, count,protocolNum},
	  protocols = GenerateProtocols[outcomes];
	  
	  count=0;
	  protocolNum=Length[protocols];
	  Monitor[Do[
	    result = CalculateMeasures[protocol, outcomes];
	    protocols[protocol] = result;
	    count+=1;
	    , {protocol, Keys[protocols]}], 
	    Row[{ProgressIndicator[count,{1,protocolNum}],StringForm["``\[LessEqual]``",count,protocolNum]}," "]];
	  Print[StringForm["Calculated `` results!",protocolNum]];
	  protocols]
CalculateAllProtocols::usage =
"
outcomes (Association[List[Number] -> Matrix]): An association where the keys are the outcomes and the values are the corresponding POVMs

Returns (Association[List[number] -> Association[...]]): An association with keys corresponding to all possible combinations of outcomes and calculated results in the values

Calculate the results for all possible protocols given an association of outcomes and their corresponding POVMs.
";
  
FindLowestMeasure[protocols_]:=Module[{lowestMeasure=Infinity, lowestMeasureProtocols},
	
	Do[
		If[protocol["measureValue"] <= lowestMeasure, 
		If[protocol["measureValue"] < lowestMeasure, Module[{}, lowestMeasureProtocols={protocol};
													 lowestMeasure=protocol["measureValue"];], 
			lowestMeasureProtocols=Append[lowestMeasureProtocols, protocol]]]
		, {protocol, protocols}];
	lowestMeasureProtocols]
FindLowestMeasure::usage = 
"
protocols (Association[List[number] -> Association[...]]): An association with keys corresponding to all possible combinations of outcomes and calculated results in the values

Returns (List[Association[...]]): A list of all the protocols with the lowest measure value

Search all calculate protocols and find the ones with the lowest figure of merit.
";

BraKetExpectation[braket_, operator_]:=ConjugateTranspose[braket] . operator . braket
BraKetExpectation::usage = 
"
braket (List[Numbers]): Vector for calculating expectation, should be of the form {a, b, ...}
operator (Matrix): Matrix to find expectation value from

Returns (Number): Expectation value
";

PartialTrace[\[Rho]_, system_, dimensions_] := Module[{operatorList, basis, operatorBasis, operator, systemDims, sum, traceOperator},
	operatorList = IdentityMatrix /@ dimensions;
	systemDims = dimensions[[system]];
	
	sum = ConstantArray[0, Dimensions[\[Rho]]/systemDims];
	Do[operatorBasis = operatorList;
		operatorBasis[[system]] = ConjugateTranspose[{UnitVector[systemDims, basis]}] (*This is given as a row vector so need to conjugate transpose*);
		traceOperator = KroneckerProduct @@ operatorBasis;
		sum += BraKetExpectation[traceOperator, \[Rho]],
	{basis, Range[dimensions[[system]]]}];
	sum]
PartialTrace::usage = 
"
\[Rho] (Matrix): Density matrix to trace over
system (Number): Identity of system to trace over. eg. trace over qubit 1 or qubit 2
dimensions (List[Number]): Dimensions of the tensor product. eg 2 qubits and a 3 level system would be {2,2,3}

Returns (Matrix): Partial trace of \[Rho]
";

EntropyOfEntanglement[state_, trace_:{1}, dimensions_:{2,2}] := Module[{\[Rho], system, \[Lambda]},
	\[Rho] = Chop[Matrixify[Normalize[state]]];
	Do[\[Rho]=Chop[PartialTrace[\[Rho], system, dimensions]], {system, trace}];
	\[Lambda] = DeleteCases[Chop[Eigenvalues[\[Rho]]], 0]; (*Only sum non-zero eigen values*)
	-Chop[Total[\[Lambda]*Log[\[Lambda]]]]]
EntropyOfEntanglement::usage = 
"
state (List[Number]): Vector of state
trace (Number): System to trace over. For a 2 qubit system it doesnt matter which

Returns (Number): Entropy of Entanglement

Entropy of entanglement between 2 systems. If 0 then the systems are fully separable.
";

PlotProtocol[protocol_, dimensions_] := Module[{grid, outcome}, 
	grid = ConstantArray[0, dimensions];
	
	Do[
		grid[[outcome[[1]],outcome[[2]]]]=1,
		{outcome, protocol}];
		MatrixPlot[grid, ColorFunction->"Monochrome"]]
PlotProtocols::usage = 
"
protocol (List[List[Number]]): A list containing the outcomes for an assignment
dimensions (List[Number]): Dimensions of the grid  

Returns (Plot): A matrix plot of the assignment outcomes. 

Outcomes in the assignment are black while outcomes not in the assignment are white.
";
	
ProbabilityGrid[state_, outcomes_] := Module[{grid, protocol},

	grid = ConstantArray[0, Dimensions[Values[outcomes][[1]]]];

	Do[
		grid[[protocol[[1]],protocol[[2]]]] = Chop[BraKetExpectation[Normalize[state], SumPovms[{protocol}, outcomes]]], {protocol, Keys[outcomes]}];
	grid
	]
ProbabilityGrid::usage = 
"
state (List[Number]): A state 
outcomes (Association[List[Number] -> Matrix]): An association where the keys are the outcomes and the values are the corresponding POVMs

Returns (Plot): A matrix plot of the probabilities

Plot the probabilities of each outcome for a given state.
";
		
ProbabilityRatioGrid[state1_, state2_, outcomes_] := Module[{grid, protocol, measureOperator, probability1, probability2},
	grid = ConstantArray[0, Dimensions[Values[outcomes][[1]]]];

	Do[
		measureOperator = SumPovms[{protocol}, outcomes];
		probability1 = Chop[BraKetExpectation[Normalize[state1], measureOperator]];
		probability2 = Chop[BraKetExpectation[Normalize[state2], measureOperator]];
		grid[[protocol[[1]],protocol[[2]]]] = Log[probability2/probability1], {protocol, Keys[outcomes]}];

	grid
	]
ProbabilityRatioGrid::usage = 
"
state1 (List[Number]): A state vector
state2 (List[Number]): A state vector
outcomes (Association[List[Number] -> Matrix]): An association where the keys are the outcomes and the values are the corresponding POVMs

Returns (Plot): A matrix plot of the probability ratios. 

Plot the ratio of probabilities for 2 different states for each outcome. Generally this will be the y/x, the probability of measuring 0 for input state |0> and the probability of measuring 0 for input state |1>.
";
		
PlotProtocolVerbose[protocol_, outcomes_, ratio_:False] := Module[{mat, protocolGrid, probabilityGrid},
(*Probability grid is calcualted assuming we are in the boring case where the |Subscript[\[Psi], 0]> = |Subscript[\[Lambda], max]> and |Subscript[\[Psi], 1]> = |Subscript[\[Lambda], min]>*)

	protocolGrid = ConstantArray[0, Dimensions[Values[outcomes][[1]]]];
	If[
		ratio, 
		probabilityGrid = ProbabilityRatioGrid[protocol["eigenVectors"][[1]], protocol["eigenVectors"][[4]], outcomes],
		probabilityGrid = ProbabilityGrid[protocol["eigenVectors"][[1]], outcomes]
	];
	probabilityGrid = (probabilityGrid-Min[probabilityGrid])/Max[probabilityGrid-Min[probabilityGrid]];
	
	Do[
		protocolGrid[[measure[[1]], measure[[2]]]]=1,
		{measure, protocol["protocol"]}];

	MatrixPlot[probabilityGrid, Mesh -> All,ColorFunction->"LightTemperatureMap",ColorFunctionScaling->False,
	      Epilog -> MapIndexed[Text[Style[Round[#, .01], Large], {#2[[2]], Length[protocolGrid]+1-#2[[1]]} - .5] &, protocolGrid, {2}]]]
PlotProtocolVerbose::usage = 
"
protocols (Association[List[number] -> Association[...]]): An association with keys corresponding to all possible combinations of outcomes and calculated results in the values
outcomes (Association[List[Number] -> Matrix]): An association where the keys are the outcomes and the values are the corresponding POVMs
ratio (Bool): True for plotting probability ratio, false for plotting absolute probability 

Returns (Plot): A matrix plot of the probability ratios or the absolute probabilities with an overlay of the assignment outcomes.

Given a protocol result, plot either the absolute probabilities or the probability ratio of y/x (depending on the ratio bool) and overlay with numbers where 1 means the outcome is in the assignment and zero means it is not.
";
	      
KeyToList[protocols_, keyFunc_]:= Module[{list},
	list = {};
	Do[list=Append[list,keyFunc[out]],{out,protocols}];
	list]
KeyToList::usage = 
"
protocols (Association[List[number] -> Association[...]]): An association with keys corresponding to all possible combinations of outcomes and calculated results in the values
keyFunc (Function): A function to be used as a key

Returns (List[...]): List of values depending on the keyFunc

Loop through all protocol results and apply the keyFunc to each one to extract values.
";
	
SplitAndSplice[state_]:=KroneckerProduct[PartialTrace[state, 2, {2,2}], PartialTrace[state, 1, {2,2}]]/Tr[state]
SplitAndSplice::usage = 
"
state (Matrix): Density matrix

Returns (Matrix): Split and splice density matrix

Expects a 2 qubit system. Partial trace out each system then tensor product back together to check if we get the same result back. This means we have a separable system.
";

NearTo[A_, B_, round_:-10]:=(Round[#1, 10^(round)] & /@ A)== (Round[#1, 10^(round)] & /@ B)
NearTo::usage = 
"
A (Number)
B (Number)
round (Number): magnitude of 10

Returns (Bool): Returns True if A and B are with round of each other. Else returns False
";

CheckSeparability[protocol_, outcomes_] := Module[{Q, Qsplice, vect, productStateTruths, \[Rho], \[Rho]splice},
	Q= Chop[N[SumPovms[protocol["protocol"], outcomes]]];
	Qsplice = SplitAndSplice[Q];
	productStateTruths = {NearTo[Q, Qsplice]};

	Do[
		\[Rho] = Chop[N[Matrixify[Normalize[vect]]]];
		\[Rho]splice = Chop[SplitAndSplice[\[Rho]]];
		AppendTo[productStateTruths, NearTo[\[Rho], \[Rho]splice, -9]], {vect, protocol["eigenVectors"]}
	];

	productStateTruths
]
CheckSeparability::usage = 
"
protocols (Association[List[number] -> Association[...]]): An association with keys corresponding to all possible combinations of outcomes and calculated results in the values
outcomes (Association[List[Number] -> Matrix]): An association where the keys are the outcomes and the values are the corresponding POVMs

Returns (List[Bool]): List of bools

Return a list of bools which state whether the POVM and eigenvectors of a protocol are separable or not.
";

CalculateProbabilities[outcomes_, state_] := Module[{probabilities},
	probabilities=<||>;
	Do[
	probabilities[key] = BraKetExpectation[state, outcomes[key]]
	, {key, Keys[outcomes]}];
	probabilities
]
CalculateProbabilities::usage = 
"
protocols (Association[List[Number] -> Association[...]]): An association with keys corresponding to all possible combinations of outcomes and calculated results in the values
state (List[Number]): State vector

Returns (Association[List[Number] -> Number]): Probabilities of a each outcome for the given input state
";

CalculateMeasuresProbability[protocol_, povmProbabilities_, \[Alpha]_, \[Beta]_] := 
 Module[{totalProb, x, y, measureValue, result},
	
  totalProb = 0;
  Do[totalProb += povmProbabilities[outcome],{outcome, protocol}];
      x = totalProb/.{\[Alpha]-> 1, \[Beta] -> 0};
      y = totalProb/.{\[Alpha]-> 0, \[Beta] -> 1};
  measureValue = MeasureAbsolute[x, y];
 
  result = <|"protocol" -> protocol, "measureValue" -> measureValue, "probability"->FullSimplify[totalProb],
    "measureFunction" -> "MeasureAbsolute", "xY" -> {x, y}|>;
  result
  ]
CalculateMeasuresProbability::usage = 
"
protocol (List[List[number]]): A list containing the outcomes for an assignment
povmProbability (Association[List[Number] -> Number]): An association where the keys are the outcomes and the values are probabilities of those outcomes (for a given state)
\[Alpha] (Variable): |0> state amplitude
\[Beta] (Variable): |1> state amplitude

Returns (Association[List[number] -> Association[...]]): An association with keys corresponding to all possible combinations of outcomes and calculated results in the values

This is like CalculateMeasures but works using probabilities for a given state and subsituted in \[Alpha] and \[Beta] to find x and y. Used to calculate if given a unitary, instead of optimising for a unitary.
";

CalculateAllProtocolsProbability[outcomes_, state_, \[Alpha]_, \[Beta]_] := 
 Module[{protocols, protocol, result, count,protocolNum, probabilities},
	  protocols = GenerateProtocols[outcomes];
      probabilities = CalculateProbabilities[outcomes, state];
	  
	  count=0;
	  protocolNum=Length[protocols];
	  Monitor[Do[
	    result = CalculateMeasuresProbability[protocol, probabilities, \[Alpha], \[Beta]];
	    protocols[protocol] = result;
	    count+=1;
	    , {protocol, Keys[protocols]}], 
	    Row[{ProgressIndicator[count,{1,protocolNum}],StringForm["``\[LessEqual]``",count,protocolNum]}," "]];
	  Print[StringForm["Calculated `` results!",protocolNum]];
	  protocols]
CalculateAllProtocolsProbability::usage = 
"
outcomes (Association[List[Number] -> Matrix]): An association where the keys are the outcomes and the values are the corresponding POVMs
state (List[Number]): A state vector
\[Alpha] (Variable): |0> state amplitude
\[Beta] (Variable): |1> state amplitude

Returns (Association[List[number] -> Association[...]]): An association with keys corresponding to all possible combinations of outcomes and calculated results in the values

This is like CalculateAllProtocols but works using probabilities. Used to calculate if given a unitary, instead of optimising for a unitary.
";
