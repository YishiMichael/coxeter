pub mod cyclotomic;
pub mod graded_algebra;
pub mod polynomial;
pub mod prufer_monoid;

#[cfg(test)]
mod tests {
    use super::graded_algebra::GradedAlgebraMeta;

    use super::cyclotomic::Cyclotomic;
    use super::cyclotomic::CyclotomicMeta;

    #[test]
    fn test_cyclotomic() {
        assert!(Cyclotomic::<i32>::root_of_unity(3, 0) == Cyclotomic::<i32>::root_of_unity(2, 0));
        dbg!(CyclotomicMeta::annihilate(
            Cyclotomic::<i32>::root_of_unity(3, 0)
                + Cyclotomic::<i32>::root_of_unity(3, 1)
                + Cyclotomic::<i32>::root_of_unity(3, 2)
        ));
        assert!(
            -Cyclotomic::<i32>::root_of_unity(3, 0)
                == Cyclotomic::<i32>::root_of_unity(3, 1) + Cyclotomic::<i32>::root_of_unity(3, 2)
        );
        assert!(
            -Cyclotomic::<i32>::root_of_unity(2, 0)
                == Cyclotomic::<i32>::root_of_unity(6, 2) + Cyclotomic::<i32>::root_of_unity(3, 2)
        );
        dbg!(CyclotomicMeta::annihilate(
            Cyclotomic::<i32>::root_of_unity(5, 1)
                + Cyclotomic::<i32>::root_of_unity(5, 2)
                + Cyclotomic::<i32>::root_of_unity(5, 3)
                + Cyclotomic::<i32>::root_of_unity(5, 4)
                - Cyclotomic::<i32>::root_of_unity(3, 1)
                - Cyclotomic::<i32>::root_of_unity(3, 2)
        ));
        assert!(
            Cyclotomic::<i32>::root_of_unity(5, 1)
                + Cyclotomic::<i32>::root_of_unity(5, 2)
                + Cyclotomic::<i32>::root_of_unity(5, 3)
                + Cyclotomic::<i32>::root_of_unity(5, 4)
                == Cyclotomic::<i32>::root_of_unity(3, 1) + Cyclotomic::<i32>::root_of_unity(3, 2)
        );
    }
}
