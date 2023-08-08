(* ::Package:: *)

(*Trine*)

phi1Trine={{1},{0}};
phi2Trine={{1/2},{Sqrt[3]/2}};
phi3Trine={{1/2},{-Sqrt[3]/2}};

(*Defining the matrices M0,M1,and M2*)
povmsTrine={2/3*Matrixify[phi1Trine],
2/3*Matrixify[phi2Trine],
2/3*Matrixify[phi3Trine]};

outcomesTrine = PovmPermutations[povmsTrine];

(*Tetra*)

phi1Tetra={1,0};
phi2Tetra={Rationalize[1/Sqrt[3]],Rationalize[Sqrt[2/3]]};
phi3Tetra={Rationalize[1/Sqrt[3]],E^(I*2*Pi/3)*Rationalize[Sqrt[2/3]]};
phi4Tetra={Rationalize[1/Sqrt[3]],E^(-I*2*Pi/3)*Rationalize[Sqrt[2/3]]};

povmsTetra = { 1/2*Matrixify[phi1Tetra],
		1/2*Matrixify[phi2Tetra],
		1/2*Matrixify[phi3Tetra],
		1/2*Matrixify[phi4Tetra]};

outcomesTetra = PovmPermutations[povmsTetra];


MeasureAbsolute[x_, y_] := ((1 - x)^2 + y^2)/(1 - x + y)

OptimalBasisValues[eigVals_] := 
 Module[{x, 
	   y, \[Lambda]min = Min[eigVals], \[Lambda]max = Max[eigVals]},
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
  
PovmPermutations[povms_] := 
 Module[{outcomes = <||>}, 
	  Do[outcomes[{i, j}] = KroneckerProduct[povms[[i]], povms[[j]]];
	   , {i, Length[povms]}, {j, Length[povms]}]; 
	  outcomes]

Matrixify[v_] := Outer[Times, v, ConjugateTranspose[v]]

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
  
SumPovms[protocol_, povmOutcomes_] := Module[{sum},
(*Use the first povm outcome assuming all povms have the same dimension*)
	  sum = ConstantArray[0, Dimensions[Values[outcomes][[1]]]];
	  
	  Do[sum += povmOutcomes[outcome]
	   , {outcome, protocol}];
	  sum]

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
  
FindLowestMeasure[protocols_]:=Module[{lowestMeasure=Infinity, lowestMeasureProtocols},
	
	Do[
		If[protocol["measureValue"] <= lowestMeasure, 
		If[protocol["measureValue"] < lowestMeasure, Module[{}, lowestMeasureProtocols={protocol};
													 lowestMeasure=protocol["measureValue"];], 
			lowestMeasureProtocols=Append[lowestMeasureProtocols, protocol]]]
		, {protocol, protocols}];
	lowestMeasureProtocols]

BraKetExpectation[braket_, operator_]:=ConjugateTranspose[braket] . operator . braket

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

EntropyOfEntanglement[state_, trace_:{1}, dimensions_:{2,2}] := Module[{\[Rho], system, \[Lambda]},
	\[Rho] = Chop[Matrixify[Normalize[state]]];
	Do[\[Rho]=Chop[PartialTrace[\[Rho], system, dimensions]], {system, trace}];
	\[Lambda] = DeleteCases[Chop[Eigenvalues[\[Rho]]], 0]; (*Only sum non-zero eigen values*)
	-Chop[Total[\[Lambda]*Log[\[Lambda]]]]]

PlotProtocol[protocol_, dimensions_] := Module[{grid, outcome}, 
	grid = ConstantArray[0, dimensions];
	
	Do[
		grid[[outcome[[1]],outcome[[2]]]]=1,
		{outcome, protocol}];
		MatrixPlot[grid, ColorFunction->"Monochrome"]]
	
ProbabilityGrid[state_, outcomes_] := Module[{grid, protocol},
	grid = ConstantArray[0, Dimensions[Values[outcomes][[1]]]];

	Do[
		grid[[protocol[[1]],protocol[[2]]]] = Chop[BraKetExpectation[Normalize[state], SumPovms[{protocol}, outcomes]]], {protocol, Keys[outcomes]}];
		grid]
		
PlotProtocolVerbose[protocol_, outcomes_] := Module[{mat, protocolGrid, probabilityGrid},

	protocolGrid = ConstantArray[0, Dimensions[Values[outcomes][[1]]]];
	probabilityGrid = ProbabilityGrid[protocol["eigenVectors"][[1]], outcomes];
	
	Do[
		protocolGrid[[measure[[1]], measure[[2]]]]=1,
		{measure, protocol["protocol"]}];

	MatrixPlot[probabilityGrid/Max[probabilityGrid], Mesh -> All,ColorFunction->"LightTemperatureMap",ColorFunctionScaling->False,
	      Epilog -> MapIndexed[Text[Style[Round[#, .01], Large], {#2[[2]], Length[protocolGrid]+1-#2[[1]]} - .5] &, protocolGrid, {2}]]]
