SELECT pet.name,owner, group_concat(owner),group_concat(remark separator '; ') 
FROM pet INNER JOIN event  
ON pet.name = event.name 
WHERE event.type = 'litter' 
GROUP by name,owner;
