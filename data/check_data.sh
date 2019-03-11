check_file_crop()
{
  echo "Checking ... $2"
  dd if="$1$2".petsc bs=1 skip=8 2> /dev/null | diff - from_spk/"$2"
}

check_files_crop()
{
  echo "Checking ... $2"
  for i in 00 01 02 03 04 05 06 07 08 09 10 11
  do 
    check_file_crop $1 "$2_$i"
  done
}

check_files()
{
  echo "Checking ... $2"
  for i in 00 01 02 03 04 05 06 07 08 09 10 11
  do
    echo "Checking ... $2_$i"
    diff "$1$2_$i.petsc" from_spk/"$2_$i"
  done
}

#check_files_crop forcing/boundary/ aice
#check_files_crop forcing/boundary/ hice
#check_files_crop forcing/boundary/ hsno
#check_files_crop forcing/boundary/ wind
#check_files_crop forcing/boundary/ swrad

#check_files forcing/domain/ Fe_dissolved
#check_files forcing/domain/ Ss
#check_files forcing/domain/ Ts
