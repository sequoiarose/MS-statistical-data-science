#View: show all the current pets and their owners in the system
SELECT * FROM current_owners_pets;

#Query 1: Find the name of all owners with more than 1 cat. (groupby & having)
SELECT O.first_name, O.last_name FROM Owners O, Pets P,  Pet_Owners PO
WHERE (PO.OID = O.OID) AND (PO.PID = P.PID) AND (P.species = 'cat')
GROUP BY O.OID HAVING count(*)>1;

#Query 2: Find the owner and vet names of owners whose pets have only had appointments with one vet. (groupby & having)
SELECT O.first_name, O.last_name, V.first_name, V.last_name from Owners O, Vets V, Pets P, Vet_Appointments VA, Pet_Appointments PA, Pet_Owners PO
WHERE (PO.OID = O.OID) AND (PO.PID = P.PID) AND (PA.PID=P.PID) AND (VA.EID = V.EID) AND (VA.AID = PA.AID)
GROUP BY O.OID HAVING count(V.EID=1);

#Query 3: Find the maximum age and species of the oldest pet (nested Query using all)
SELECT P.species, P.age FROM Pets P
WHERE P.age>=all (SELECT P2.age FROM Pets P2);

#Query 4: Find all pets with no vaccines (nested query using in)
SELECT * FROM Pets P
WHERE P.PID NOT IN (SELECT PV.PID FROM Pet_Vaccines PV);
