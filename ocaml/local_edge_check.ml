(* local_edge_check_diff.ml *)

open Yojson.Safe
open Yojson.Safe.Util

exception Check_error of string

let failf fmt = Printf.ksprintf (fun s -> raise (Check_error s)) fmt

let int_array_of_json_list j =
  j |> to_list |> List.map to_int |> Array.of_list

let strictly_sorted a =
  let rec aux i =
    i + 1 >= Array.length a || (a.(i) < a.(i + 1) && aux (i + 1))
  in
  aux 0

let contains_sorted a x =
  let rec aux lo hi =
    if lo >= hi then false
    else
      let mid = (lo + hi) / 2 in
      if a.(mid) = x then true
      else if a.(mid) < x then aux (mid + 1) hi
      else aux lo mid
  in
  aux 0 (Array.length a)

let sorted_difference a b =
  (* Return a \ b, assuming both arrays are sorted. *)
  let na = Array.length a in
  let nb = Array.length b in
  let out = ref [] in
  let i = ref 0 in
  let j = ref 0 in

  while !i < na && !j < nb do
    if a.(!i) = b.(!j) then begin
      incr i;
      incr j
    end else if a.(!i) < b.(!j) then begin
      out := a.(!i) :: !out;
      incr i
    end else
      incr j
  done;

  while !i < na do
    out := a.(!i) :: !out;
    incr i
  done;

  Array.of_list (List.rev !out)

let sorted_intersects a b =
  let na = Array.length a in
  let nb = Array.length b in
  let i = ref 0 in
  let j = ref 0 in
  let found = ref false in

  while (not !found) && !i < na && !j < nb do
    if a.(!i) = b.(!j) then
      found := true
    else if a.(!i) < b.(!j) then
      incr i
    else
      incr j
  done;

  !found

let parse_certificate path =
  let json = Yojson.Safe.from_file path in

  let items_json = json |> member "items" |> to_list in
  let incidents =
    items_json
    |> List.map (fun item -> item |> member "incident" |> int_array_of_json_list)
    |> Array.of_list
  in

  let neighbors =
    json |> member "neighbors" |> to_list
    |> List.map int_array_of_json_list
    |> Array.of_list
  in

  incidents, neighbors

let check_neighbors incidents neighbors =
  let n = Array.length incidents in

  if Array.length neighbors <> n then
    failf "neighbors has length %d, but items has length %d"
      (Array.length neighbors) n;

  Array.iteri
    (fun v inc ->
      if not (strictly_sorted inc) then
        failf "incident list of item %d is not strictly sorted" v)
    incidents;

  Array.iteri
    (fun v ns ->
      if not (strictly_sorted ns) then
        failf "neighbors[%d] is not strictly sorted" v;

      Array.iter
        (fun w ->
          if w < 0 || w >= n then
            failf "neighbors[%d] contains out-of-range item %d" v w;
          if w = v then
            failf "neighbors[%d] contains itself" v;

          if not (contains_sorted neighbors.(w) v) then
            failf "neighbor relation is not symmetric: %d lists %d, but not conversely" v w)
        ns)
    neighbors

let check_local_edge_test incidents neighbors =
  let n = Array.length incidents in

  for v = 0 to n - 1 do
    let iv = incidents.(v) in
    let nv = neighbors.(v) in

    (* Transposed formulation:
       I(v) ∩ I(w) ⊄ I(u) iff I(w) ∩ (I(v) \ I(u)) ≠ ∅. *)
    Array.iter
      (fun u ->
        let iu = incidents.(u) in
        let diff_vu = sorted_difference iv iu in

        Array.iter
          (fun w ->
            if w <> u then begin
              let iw = incidents.(w) in
              if not (sorted_intersects iw diff_vu) then
                failf
                  "local edge test failed at item %d: neighbor %d is not separated from neighbor %d; I(%d) ∩ I(%d) is contained in I(%d)"
                  v w u v w u
            end)
          nv)
      nv
  done

let () =
  if Array.length Sys.argv <> 2 then begin
    Printf.eprintf "usage: %s certificate.json\n" Sys.argv.(0);
    exit 2
  end;

  let path = Sys.argv.(1) in

  try
    let incidents, neighbors = parse_certificate path in
    check_neighbors incidents neighbors;

    let t0 = Unix.gettimeofday () in
    check_local_edge_test incidents neighbors;
    let t1 = Unix.gettimeofday () in

    Printf.printf "OK: local edge test passed\n";
    Printf.printf "local edge test time: %.6f seconds\n" (t1 -. t0)
  with
  | Check_error msg ->
      Printf.eprintf "ERROR: %s\n" msg;
      exit 1
  | Yojson.Json_error msg ->
      Printf.eprintf "JSON ERROR: %s\n" msg;
      exit 1
  | Type_error (msg, _) ->
      Printf.eprintf "JSON TYPE ERROR: %s\n" msg;
      exit 1
