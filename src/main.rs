

fn pbwt(n_haplotypes: &usize, n_variants: &usize, data: &Vec<usize>, min_length: &usize) -> Vec<[usize; 3]> {
    let mut ppa: Vec<usize> = (0..*n_haplotypes).collect();
    let mut div: Vec<usize> = vec![0; *n_haplotypes];
    let mut results: Vec<[usize; 3]> = Vec::new();

    for i in 0..*n_variants {
        let mut ppa_0 = Vec::<usize>::with_capacity(data.len());
        let mut ppa_1 = Vec::<usize>::with_capacity(data.len());
        let mut div_0 = Vec::<usize>::with_capacity(data.len());
        let mut div_1 = Vec::<usize>::with_capacity(data.len());
        let (mut p, mut q) = (i + 1, i + 1);
        let mut mat_0 = Vec::<usize>::with_capacity(data.len());
        let mut mat_1 = Vec::<usize>::with_capacity(data.len());

        for j in 0..ppa.len() {
            let match_start = div[j];

            if match_start + min_length > i {
                if mat_0.len() > 0 && mat_1.len() > 0 {
                    gather_results(&mut results, &i, &mat_0, &mat_1);
                }
                mat_0.clear();
                mat_1.clear();
            }

            let index_ppa = ppa[j];
            let allele = data[index_ppa * n_variants + i];

            if match_start > p {
                p = match_start;
            }
            if match_start > q {
                q = match_start;
            }

            if allele == 0 {
                ppa_0.push(index_ppa);
                div_0.push(p);
                p = 0;
                mat_0.push(index_ppa);
            } else {
                ppa_1.push(index_ppa);
                div_1.push(q);
                q = 0;
                mat_1.push(index_ppa);
            }
        }

        if mat_0.len() > 0 && mat_1.len() > 0 {
            gather_results(&mut results, &i, &mat_0, &mat_1);
        }

        ppa.clear();
        ppa.append(&mut ppa_0);
        ppa.append(&mut ppa_1);

        div.clear();
        div.append(&mut div_0);
        div.append(&mut div_1);
    }

    {
        let mut temp_data = Vec::<usize>::new();

        for index in ppa.iter() {
            temp_data.extend_from_slice(&data[(index * n_variants)..(index * n_variants + n_variants)]);
        }
        print_data(&n_haplotypes, &n_variants, &temp_data);
    }

    println!("{:?}", div);
    return results;
}

fn gather_results(results: &mut Vec<[usize; 3]>, index: &usize, mat_0: &Vec<usize>, mat_1: &Vec<usize>) {
    for hap_1 in mat_0.iter() {
        for hap_2 in mat_1.iter() {
            results.push([*index, *hap_1, *hap_2]);
        }
    }
}

fn print_data(n_haplotypes: &usize, n_variants: &usize, data: &Vec<usize>) {
    for i in 0..*n_haplotypes {
        println!("{:?}", &data[(i * n_variants)..(i * n_variants + n_variants)]);
    }
    println!("");
}

fn counting_sort(pbwt_results: Vec<[usize; 3]>, n_haplotypes: &usize) -> Result<Vec<[usize; 3]>, String> {
    let mut count_1_results: Vec<[usize; 3]>;

    if pbwt_results.len() > 0 {
        count_1_results = counting_sub_sort(pbwt_results, &n_haplotypes, 1);
    } else {
        return Err(String::from("No PBWT results."));
    }

    let mut start = 0;
    let mut current: Option<usize>;
    
    if count_1_results.len() > 0 {
        current = Some(count_1_results[0][1]);
    } else {
        return Err(String::from("No counting sort results."));
    }

    let mut count_2_slice: Vec<[usize; 3]>;
    let mut count_2_slice_results: Vec<[usize; 3]>;

    if current.is_some() {
        for i in 0..count_1_results.len() {
            if current.unwrap() != count_1_results[i][1] {
                count_2_slice = count_1_results[start..i].to_vec();

                if count_2_slice.len() > 1 {
                    count_2_slice_results = counting_sub_sort(count_2_slice, &n_haplotypes, 2);
                    count_1_results.splice(start..i, count_2_slice_results);
                }

                start = i;
                current = Some(count_1_results[i][1]);
            }
        }

        if start != count_1_results.len() - 1 {
            let end_index = count_1_results.len();
            count_2_slice = count_1_results[start..end_index].to_vec();

            if count_2_slice.len() > 1 {
                count_2_slice_results = counting_sub_sort(count_2_slice, &n_haplotypes, 2);
                count_1_results.splice(start..end_index, count_2_slice_results);
            }
        }
    } else {
        return Err(format!("Invalid current value: {}", current.unwrap()));
    }
    Ok(count_1_results)
}

fn counting_sub_sort(pbwt_results: Vec<[usize; 3]>, n_haplotypes: &usize, pair_index: usize) -> Vec<[usize; 3]> {
    let mut count_vec: Vec<usize> = vec![0; *n_haplotypes];
    let mut count_results: Vec<[usize; 3]> = vec![[0, 0, 0]; pbwt_results.len()];

    for result in pbwt_results.iter() {
        count_vec[result[pair_index]] += 1;
    }

    for i in 1..count_vec.len() {
        count_vec[i] += count_vec[i - 1];
    }

    for result in pbwt_results {
        count_results[count_vec[result[pair_index]] - 1] = result;
        count_vec[result[pair_index]] -= 1;
    }

    return count_results
}

fn main() {
    const N_HAPLOTYPES: usize = 8;
    const N_VARIANTS: usize = 6;
    let data = vec![0, 1, 0, 1, 0, 1,
                    1, 1, 0, 0, 0, 1,
                    1, 1, 1, 1, 1, 1,
                    0, 1, 1, 1, 1, 0,
                    0, 0, 0, 0, 0, 0,
                    1, 0, 0, 0, 1, 0,
                    1, 1, 0, 0, 0, 1,
                    0, 1, 0, 1, 1, 0];
    const MIN_LENGTH: usize = 3;

    if N_HAPLOTYPES * N_VARIANTS != data.len() {
        panic!("Dimensions mismatch: haplotypes: {}, variants: {}, data: {}", N_HAPLOTYPES, N_VARIANTS, data.len());
    }

    print_data(&N_HAPLOTYPES, &N_VARIANTS, &data);
    let mut results = pbwt(&N_HAPLOTYPES, &N_VARIANTS, &data, &MIN_LENGTH);
    println!("{:?}", results);
    results = match counting_sort(results, &N_HAPLOTYPES) {
        Ok(ok) => ok,
        Err(err) => panic!("{:?}", err)
    };
    println!("{:?}", results);
}

#[cfg(test)]
mod test {
    #[test]
    fn test_pbwt_integration() {
        use main;

        main();
        assert!(true);
    }

    #[test]
    fn test_pbwt() {
        use pbwt;

        const N_HAPLOTYPES: usize = 8;
        const N_VARIANTS: usize = 6;
        const MIN_LENGTH: usize = 3;
        let data = vec![0, 1, 0, 1, 0, 1,
                        1, 1, 0, 0, 0, 1,
                        1, 1, 1, 1, 1, 1,
                        0, 1, 1, 1, 1, 0,
                        0, 0, 0, 0, 0, 0,
                        1, 0, 0, 0, 1, 0,
                        1, 1, 0, 0, 0, 1,
                        0, 1, 0, 1, 1, 0];

        assert_eq!(data.len(), N_HAPLOTYPES * N_VARIANTS);

        let results = pbwt(&N_HAPLOTYPES, &N_VARIANTS, &data, &MIN_LENGTH);
        let target = vec![[4, 4, 5], [4, 0, 7], [5, 4, 1], [5, 4, 6], [5, 3, 2]];

        assert_eq!(target, results);
    }

    #[test]
    fn test_counting_sort() {
        use counting_sort;

        const N_HAPLOTYPES: usize = 8;
        let data = vec![[4, 4, 5], [4, 0, 7], [5, 4, 1], [5, 4, 6], [5, 3, 2]];
        let results = match counting_sort(data, &N_HAPLOTYPES) {
            Ok(ok) => ok,
            Err(err) => panic!("{:?}", err)
        };
        let target = vec![[4, 0, 7], [5, 3, 2], [5, 4, 1], [4, 4, 5], [5, 4, 6]];

        assert_eq!(target, results);
    }
}