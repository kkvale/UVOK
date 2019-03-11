check_file()
{
  echo "Checking ... $1"
  dd if=forcing/boundary/"$1".petsc bs=1 skip=8 2> /dev/null | diff - from_spk/"$1"
}

for i in 00 01 02 03 04 05 06 07 08 09 10 11
do 
  check_file "aice_$i"
done

